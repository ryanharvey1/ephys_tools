function CheckEegStates(FileBase,varargin)
% function CheckEegStates(FileBase,State, AuxData,FreqRange,Channels,Window, Action, Overwrite)
% starts the gui for browsing through spectral representation of the eeg and
% segmenting into states .. has lot's of features - try it. to speed up
% precompute the spectrograms using EegSegmentation function. 
% if no Channels are specified, and you didn't create FileBase.eegseg.par
% file with channels (starting from 0) that you want to use, it will start
% a dialog, which runs away under new version of matlab ..don't know why.
% Window - in sec size of the spectral window
% FreqRange =[Fmin Fmax] - range of freq over which the spectrum will be
% computed
% AuxData may contain the additional data to plot
% AuxData has to be cell array where each row is : xaxis, yaxis, data, display_func
% obviously the xaxis (time) has to match the full time of the file, and be in the
% same units as spectrogram - seconds
% so far display_func is 'plot' and 'imagesc'
% in case of 'plot' the yaxis can be empty. 
% Action : 'display' or 'compute'
% The following keys switch (on/off) the following behavior :
%     t (default) - to move the current position pointer and display the corresponding trace of eeg', ...
%     n - add new region first click set''s one border, second - the other. new region is displayed as green/red lines ..etc',...
%     m - move one border. click next to the border and drag it where you want. ',..
%     toggles t,n,m are not really neded if you have 3-button mouse, much
%     faster is to use the mouse buttons as described below.
%     d - delete the border. Just click on the border
%     z - toggles zoom state. Mouse :left/right btn - zoom in/out the position of the mouse. Keybd: f - resets x axis to max. arrows up/down- zoom in/out the current position',...
%     c - and then  up/down - rescale the color axis in spectrograms',...
%     f - and then  up/down - zoom in/out the y axis in spectrograms',...
%     w - and then up/down - decrease/increas the size of the window for the traces display ',...
% These keys do
%     h - shows help box
%     u - update periods. Removes deleted and reorders all lines, pairs consequituive ones. Colorizes blue -beg, red -end',...
%     s - save the results in the file (currently ASCII two columnt, 1250 s.r.', ...
%     l - load segmentation periods from a file ', ...
%     arrows left/right - move in the spectrogram', ...
%     space - move to the closest fromt the right beginning of a period', ....
%     a - automatics segementation tool. To use it mark boundaries of a
%     segment you want to process automatically, click inside and press
%     'a'. A menu to select passband/stopband and channel to use will show
%     up. The stopband can be more than one set of bands. The algorythm
%     will compute ratio of power in pass/stop bands and perform HMM two
%     stage fit. New segments that correspond to the state with high values
%     of the power ratio will appear. You may have to wait a bit.
% Mouse click in trace (default) mode does: 
%     left - view trace at the point of click
%     right - add new border, middle - move the border


[State,AuxData,FreqRange,Channels,Window,Action, Overwrite] = ...
    DefaultArgs(varargin,{[],[],[1 100],[],1,'display',0});

% auxil. struct for gui
global gCheckEegStates
gCheckEegStates = struct;

%constants
UseRulesBefore = 0; % flag to switch heuristics for periods cleaning after automatic segmentation
MinLen=5; %seconds

Par = LoadXml([FileBase '.xml']);

if isfield(Par,'lfpSampleRate')
    eSampleRate = Par.lfpSampleRate;
else
    eSampleRate = 1250;
end

if ~isempty(State)
    % load segmentation results
    if exist([FileBase '.sts.' State],'file')
        Per = load([FileBase '.sts.' State]);
%     elseif FileExists([FileBase '.states.res'])
%         Per = SelectStates(FileBase,State,eSampleRate*MinLen);
    else
        Per = [];
    end       
    if isempty(Per)
        fprintf('No segments detected in %s state\n',State)
        %return;
        Per = [];
    end
else
    Per = [];
    State = '';
end

if isempty(Channels) && exist([FileBase '.eegseg.par'],'file');
    Channels = load([FileBase '.eegseg.par'])+1;
    %         Channels = GetChannels(FileBase, {{'h','c'}});
end
if ~isempty(Channels) && ~ exist([FileBase '.eegseg.par'],'file');
    fprintf('You selection of Channels is stored in the file %s for future usage, so that next time you don''t have to pass Channels\n',...
        [FileBase '.eegseg.par']);
    msave([FileBase '.eegseg.par'],Channels-1);
end
Compute=0;
if exist([FileBase '.eegseg.mat'],'file') 
    switch Action
        case 'display'
            if Overwrite ==0
                over='No';
            else
                over = questdlg('Do you want to overwrite existing spectrogram?','Overwrite');
            end
            switch over
                case 'No'
                    load([FileBase '.eegseg.mat']); % load whitened spectrograms from EegSegmentation results
                case 'Yes'
                    Compute=1;
            end
        case 'compute'
            if Overwrite
                Compute=1;
            else
                return;
            end
    end
end
%keyboard
if ~exist([FileBase '.eegseg.mat'],'file') || Compute==1

    if isempty(Channels)
        ch = inputdlg({'Enter channels to use'},'Channels',1,{'1'});
        Channels = str2num(ch{1});
    end
 
    % now compute the spectrogram
    if exist([FileBase '.eeg'],'file')
        try 
            Eeg = LoadBinary([FileBase '.eeg'],'channels',Channels,'nChannels',str2num(Par.nChannels))';
        catch
            Eeg = LoadBinary([FileBase '.eeg'],'channels',Channels,'nChannels',Par.nChannels)';
        end
    elseif exist([FileBase '.eeg.0'],'file')
        Eeg = bload([FileBase '.eeg.0'],[1 inf]);
        
    else
        error('no eeg file or eeg.0 file! \n');
    end
    fprintf('computing spectrograms, may take time ... \n');
    %nFFT = 2^round(log2(2^11)); %compute nFFT according to different sampling rates
    SpecWindow = 2^round(log2(Window*eSampleRate));% choose window length as power of two
    nFFT = SpecWindow*4;
    weeg = WhitenSignal(Eeg,eSampleRate*2000,1);
    [y,f,t]=mtcsglong(weeg,nFFT,eSampleRate,SpecWindow,[],2,'linear',[],FreqRange);
    save([FileBase '.eegseg.mat'],'y','f','t','Channels','-v6');
end

if strcmp(Action,'compute')
    return;
end
t = (t(2)-t(1))/2 +t;
   
% computer the/del ratio and detect transitions automatically - not used at
% the momnet, maybe later
%[thratio] = TDRatioAuto(y,f,t,MinLen);
%[thratio, ThePeriods] = TDRatioAuto(y,f,t,MinLen);

%now apply the rules to filter out junk states or make continuous periods
% to be implemented later
if UseRulesBefore
    switch State
        case 'REM'

    end
end

% fill the global structure for future use

if ~exist([FileBase '.eeg'],'file') && FileExists([FileBase '.eeg.0'])
    gCheckEegStates.EegFile  ='eeg.0';
    gCheckEegStates.Channels =1;
    gCheckEegStates.nChannels = 1;
    gCheckEegStates.nSamples = FileLength([FileBase '.eeg.0'])/2;
else
    gCheckEegStates.EegFile  ='eeg';
    gCheckEegStates.Channels = Channels;
    gCheckEegStates.nChannels = length(Channels);
    gCheckEegStates.nSamples = FileLength([FileBase '.eeg'])/Par.nChannels/2;
end

nAuxData = max(size(AuxData,1));

gCheckEegStates.FileBase = FileBase;
gCheckEegStates.Par = Par;
gCheckEegStates.State = State;
gCheckEegStates.State2 = ''; %State2, Periods2, lh2 are used for a secondary (not editable) state
gCheckEegStates.t = 10; %is seconds
gCheckEegStates.eFs = eSampleRate;
gCheckEegStates.trange = [t(1) t(end)];
gCheckEegStates.Periods = Per/eSampleRate; % in seconds
gCheckEegStates.Periods2 = []; % in seconds
gCheckEegStates.Mode = 't';
gCheckEegStates.nPlots=gCheckEegStates.nChannels+1+nAuxData;
gCheckEegStates.lh=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.lh2=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.Window = Window*eSampleRate*2;
gCheckEegStates.SelLine=[];
gCheckEegStates.cposh=cell(gCheckEegStates.nPlots,1);
gCheckEegStates.FreqRange = [min(f) max(f)];
gCheckEegStates.newl=[];
gCheckEegStates.tstep = t(2)-t(1);
gCheckEegStates.coolln = [];
gCheckEegStates.LastBut = 'normal';
gCheckEegStates.nAuxData = nAuxData;
gCheckEegStates.nAuxData = nAuxData;
gCheckEegStates.RefPoints = [];
gCheckEegStates.Spec = y;
gCheckEegStates.Freq = f;
gCheckEegStates.Times = t;
gCheckEegStates.SubPlots = zeros(gCheckEegStates.nPlots,1);
gCheckEegStates.ThetaRatioPlot = [];
Data=AuxData(1,:);
gCheckEegStates.AuxTimes = Data{1};
gCheckEegStates.AuxData = Data{3};

if nAuxData>0
    gCheckEegStates.AuxDataType = AuxData(:,4);
end
% create and configure the figure
gCheckEegStates.figh = figure('ToolBar','none');
%set(gCheckEegStates.figh, 'Position', [3 828 1276 620]); %change Postion of figure if you have low resolution
set(gCheckEegStates.figh, 'Name', 'CheckEegStates : traces');
set(gCheckEegStates.figh, 'NumberTitle', 'off');


% put the uitoolbar and uimenu definitions here .. may require rewriting
% some callbacks as actions rather then cases of actions (e.g. key
% pressing)


%now do the plots

for ii=1:gCheckEegStates.nChannels
    gCheckEegStates.SubPlots(ii) = subplot(gCheckEegStates.nPlots,1,ii);
    imagesc(t,f,log(sqrt(y(:,:,ii)))');axis xy; ylim([max(0,FreqRange(1)) min(FreqRange(2),20)]);
    hold on
    if ii==1
        title('Spectrogram'); ylabel('Frequency (Hz)');
    end
end

if nAuxData>0
    for ii=[1:nAuxData]
        gCheckEegStates.SubPlots(ii+gCheckEegStates.nChannels) = subplot(gCheckEegStates.nPlots,1,ii+gCheckEegStates.nChannels);
        DisplayAuxData(AuxData(ii,:));
        xlim(gCheckEegStates.trange);
        hold on
        
    end
end


%  subplot(gCheckEegStates.nPlots,1,2)
%  plot(t,thratio);axis tight;
%  set(gca,'YTick',[]);
%  hold on
%  ylabel('Theta/Delta raio'); xlabel('Seconds');

subplot(gCheckEegStates.nPlots,1,gCheckEegStates.nPlots)
CheckEegStates_aux('traces'); % plot the current eeg traces
hold on

%plots lines
if ~isempty(Per)
CheckEegStates_aux('lines');
end
% assign functions for mouse and keyboard click
set(gCheckEegStates.figh,'WindowButtonDownFcn','CheckEegStates_aux(''mouseclick'')');
set(gCheckEegStates.figh,'KeyPressFcn', 'CheckEegStates_aux(''keyboard'')');

%msave([FileBase '.states.' State],round(ThePeriods*eSampleRate));
%clear y
return




function [thratio, ThePeriods] = TDRatioAuto(y,f,t,MinLen)

%automatic theta periods detection just using the thetaratio
thfin = find(f>6 & f<9);
thfout = find(f<5 | (f> 12& f<25));
thratio = log(mean(sq(y(:,thfin,1)),2))-log(mean(sq(y(:,thfout,1)),2));

if nargout>1
    nStates =2;
    % fit gaussian mixture and to HMM - experimental version .. uses only thetaratio
    [TheState thhmm thdec] = gausshmm(thratio,nStates,1,0);

    for i=1:nStates
        thratio_st(i) = mean(thratio(TheState==i));
    end

    [dummy TheInd] = max(thratio_st);
    InTh = (TheState==TheInd);
    DeltaT = t(2)-t(1);
    MinSeg = round(MinLen/DeltaT);

    TransTime = ThreshCross(InTh,0.5,MinSeg);
    ThePeriods = t(TransTime);
end
return

function DisplayAuxData(Data)

        nEl =  size(Data,2);
        if nEl<4 
            err=1;
        elseif ~isstr(Data{4})
            err=1;
        else
            err=0;
        end
        if err 
            warning('AuxData has to be cell array where each row is : xaxis, yaxis, data, display_func');
            close 
            return;
        end
        
        switch Data{4} %switch by functions
            case 'plot'
                plot(Data{1}, Data{3});
            case 'imagesc'
                if length(Data{1})~=size(Data{3},1) & length(Data{1})~=size(Data{3},2)
                    Data{3}=Data{3}';
                end
                imagesc(Data{1},Data{2}, Data{3}');
                axis xy
            otherwise 
                error('wrong data display function');
        end
        


return