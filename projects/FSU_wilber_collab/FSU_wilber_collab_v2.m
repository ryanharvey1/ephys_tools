function data=FSU_wilber_collab_v2(path)
% FSU_wilber_collab main function for linear track anaysis
%
%   Input:
%           path: path to session data*
%   Output:
%           data: data structure containing tracking/spike data,
%                   measure results, and more
%
% *Example:
%   main folder: Alz39_31917
%   contents:   .t spike files
%               .nvt video file
%               .txt event file
%
% dependencies: ALLFUNCS.m, Mat2NLxVT.mexw64(or mac equivalent), ts.m, Data.m
%
% Ryan Harvey 2019


%% move to session directory
cd(path)

%% set up id & data struct
data=ALLFUNCS.datastruct(path);
clear path

%% check for notes
if exist('notes.mat','file')
    load('notes.mat','notes')
    data.notes=notes;
    clear notes
end
%% task parameters
data.mazetype='LinearTrack';
data.mazesize=120;

%% read & parse events
data.events=ALLFUNCS.parsewilberevents;

%% read trial ts
trial_ts=load('trial_ts.mat');
name=fieldnames(trial_ts);
data.trial_ts=trial_ts.(name{1});
clear trial_ts name

%% load lfp
[data]=ALLFUNCS.handle_LFP(data);

%% extract movement info
[ts,x,y,angles] = Nlx2MatVT([data.session_path,filesep,'VT1.nvt'],[1,1,1,1,0,0],0,1);

%% check for manual coordinate correction
if exist([data.session_path,filesep,'restrictxy.mat'],'file')
    % run the following line if tracker errors from unplugs are present
    %     ALLFUNCS.manual_trackerjumps(ts,x,y,data.events(1),data.events(2),data.session_path);
    load([data.session_path,filesep,'restrictxy.mat'],'in')
    x(in==0)=NaN;
    y(in==0)=NaN;
    clear in
end

%% calculate video sample rate
data.samplerate=ALLFUNCS.samplerate(ts);

%% fix non-detects and smooth
[xtemp,ytemp]=ALLFUNCS.FixPos(x',y',ts',round(0.1667*data.samplerate));

tempangle=wrapTo360(ALLFUNCS.fixNLXangle(angles',round(0.1667*data.samplerate)))'-90;
tempangle(tempangle<0)=tempangle(tempangle<0)+360;

data.frames=[ts',rescale(xtemp,0,data.mazesize),rescale(ytemp,0,14),tempangle];

clear xtemp ytemp ts tempangle x y angles

%% load spikes
files=dir('*.t');
data.cellid={files.name};
data.Spikes=ALLFUNCS.LoadSpikes(data.cellid);
clear files

%% fix time: normalizes time to the first tracker frame and converts ts to seconds
data=ALLFUNCS.FixTime(data);

%% limit frames to session
if ~isempty(data.events)
    data.frames=data.frames(data.frames(:,1)>=data.events(1) & data.frames(:,1)<=data.events(2),:);
end

%% loop through each cell
for i=1:length(data.Spikes)
    %% split trials up
    [data.ratemaps{i},data.occ{i},data.frames_w_spk{i},data.trial_frames_w_spk{i}]=ALLFUNCS.split_trials(data,data.Spikes{i});
    %% split up right and left
    if isfield(data,'notes')
        if contains(data.notes.trials,'rightleft')
            [data.right{i},data.left{i}]=ALLFUNCS.split_trials_RL(data.ratemaps{i},data.occ{i},data.trial_frames_w_spk{i});
        end
    end
    %% info content per trial
%     data.trials.trial_infocontent{i}=ALLFUNCS.infocontent(data.ratemaps{i},data.occ{i});
    
    %% THETA MODULATION
    [data.results.thetaindex(i),data.results.peaktheta(i),data.autocor(i,:)]=...
        ALLFUNCS.thetamodulation(data.frames_w_spk{i}(data.frames_w_spk{i}(:,end)==1,1));
end
% save data
cd ..
cd ..
save(['ProcessedData\',data.rat,'_',data.sessionID],'-struct','data','-v7.3')


%% plot session figures
close all
wilberFigures_v2(data);

figuresize=[2 42 838 924];
if isfield(data,'notes')
    if contains(data.notes.trials,'rightleft')
        figuresize=[1 41 1680 933];
    end
end
FolderName='D:\Projects\FSU_wilber_collab\Figures';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');  
  set(FigHandle,'Position',figuresize)
  print(FigHandle,'-dpng', '-r300', [FolderName,filesep,erase(FigName,{'.t',':',' '}),'.png'])
  saveas(FigHandle, [FolderName,filesep,erase(FigName,{'.t',':',' '}),'.fig'])
end
close all

end