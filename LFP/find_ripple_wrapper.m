function ripple_info = find_ripple_wrapper(data,varargin)
%find_ripple_wrapper: wrapper for findRipples from TSToolbox_Utils
%
%   Input:
%           data: ephys_tools data structure. (if you don't have this, you 
%                                               can use [] for the first input)
%           varargin:
%               figs: debugging figs, 0 or 1
%               lfp: channels x time 
%               lfp_ts: time stamps from lfp in seconds
%               speed: speed of animal (cm/sec)
%               mov_ts: time stamps from movement in seconds
%
%   Output:
%           ripple_info:
%                ripple_info.ripple_start: start time of ripple
%                ripple_info.ripple_peak: peak time of ripple
%                ripple_info.ripple_end: end time of ripple
%                ripple_info.peakNormalizedPower: power of ripple
%                ripple_info.frequency: frequency of ripple
%                ripple_info.amplitude: peak amplitude of ripple
%                ripple_info.ripple_dur: duration of ripple in seconds
%                ripple_info.unfiltered_ripple: unfiltered ripple
%                ripple_info.filtered_ripple: filtered ripple
%                ripple_info.time: time stamps during ripple
%
% dependent on buzcode.  ephys_tools\external_packages\buzcode
%
% Ryan Harvey 2019

if isstruct(data)
    % unpack input
    p = inputParser;
    p.addParameter('figs',0);
    p.addParameter('lfp',data.lfp.signal);
    p.addParameter('lfp_ts',data.lfp.ts);
    p.addParameter('frequency',data.lfp.lfpsamplerate);
    p.addParameter('speed',data.frames(:,5));
    p.addParameter('mov_ts',data.frames(:,1));
    
    p.parse(varargin{:});
    figs = p.Results.figs;
    lfp = p.Results.lfp;
    lfp_ts = p.Results.lfp_ts;
    frequency = p.Results.frequency;
    speed = p.Results.speed;
    mov_ts = p.Results.mov_ts;
else
    p = inputParser;
    p.addParameter('figs',0);
    p.addParameter('lfp',[]);
    p.addParameter('lfp_ts',[]);
    p.addParameter('frequency',[]);
    p.addParameter('speed',[]);
    p.addParameter('mov_ts',[]);
    
    p.parse(varargin{:});
    figs = p.Results.figs;
    lfp = p.Results.lfp;
    lfp_ts = p.Results.lfp_ts;
    frequency = p.Results.frequency;
    speed = p.Results.speed;
    mov_ts = p.Results.mov_ts;
end


% check for EMG channel
parts = strsplit(data.session_path,filesep);
basepath = fullfile(parts{1:find(parts == "Projects",1,'last') + 1});
if exist(fullfile(basepath,'EMG_from_LFP'),'dir')
    processedpath=strsplit(data.session_path,filesep);
    processedpath(end-2:end)=[];
    emg_file = fullfile(strjoin(processedpath,filesep),'EMG_from_LFP',...
        [data.rat,'_',data.sessionID,'_emg.mat']);
    if exist(emg_file,'file')
       emg = load(emg_file); 
       emg = emg.data;
    end
else
    emg = zeros(size(lfp,2),1);
end

for ch =1:size(lfp,1)
    signal = lfp(ch,:);
    
    % filter signal to ripple range [150 250 hz]
    signal_filtered = BandpassFilter(signal, 1000, [150 250]);
    
    [ripples] = bz_FindRipples_ephys_tools(signal',lfp_ts',...
        'EMGfilename',emg_file,...
        'passband',[150 250],...
        'frequency',frequency,...
        'thresholds',[3 7],...
        'durations',[30 350],...
        'minDuration',25,...
        'EMGThresh',0.9);
    [maps,ripple_data,stats] = bz_RippleStats(signal_filtered',lfp_ts',ripples,'frequency',frequency);
    
    bz_PlotRippleStats(maps,ripple_data,stats,'frequency',frequency);
    % locate potential ripples
%     [ripples,~] = findRipples(signal_filtered','frequency',frequency,'noise',emg);
    
    if isempty(ripples)
        continue
    end
    
    % align to session
    ripples(:,1:3) = ripples(:,1:3) - lfp_ts(1);
    
    % exclude movement
    temp_speed = interp1(mov_ts,speed,ripples(:,1:3));
    ripples = ripples(any(temp_speed < 4 | isnan(temp_speed),2),:);
    
    % max amplitude
    for r = 1:size(ripples)
        ripples(r,6) = max(signal_filtered(lfp_ts >= ripples(r,1) & lfp_ts <= ripples(r,3)));
    end
    
    % number and mean power
    n_ripples(ch) = (size(ripples,1));
    mean_power(ch) = nanmean(ripples(:,4));
    
    % store ripple info
    all_ripples{ch} = ripples;
end

% if no ripple is found
if isempty(ripples)
    ripple_info.ripple_start = NaN;
    ripple_info.ripple_peak = NaN;
    ripple_info.ripple_end = NaN;
    ripple_info.peakNormalizedPower = NaN;
    ripple_info.frequency = NaN;
    ripple_info.amplitude = NaN;
    ripple_info.ripple_dur = NaN;
    ripple_info.unfiltered_ripple = NaN;
    ripple_info.filtered_ripple = NaN;
    ripple_info.time = NaN;
    return
end
clear ripples
for ch = 1:length(all_ripples)
    % re-bandpass filter
    signal_filtered = BandpassFilter(lfp(ch,:), 1000, [150 250]);
    
    ripples= all_ripples{ch};
    
    % pull out unfiltered and filtered ripples from lfp
    for r = 1:size(ripples,1)
        idx = lfp_ts >= ripples(r,1) & lfp_ts <= ripples(r,3);
        unfiltered_ripple{ch}{r} = signal(idx);
        filtered_ripple{ch}{r} = signal_filtered(idx);
        time{ch}{r} = lfp_ts(idx);
    end
    % get ripple duration in sec
    ripple_dur{ch} = ripples(:,3) - ripples(:,1);
end

% % locate channel with highest mean power and number of ripples
% [~,ch] = max(mean_power.*n_ripples);
% ripples = all_ripples{ch};
% signal = lfp(ch,:);
% 
% % re-bandpass filter
% signal_filtered = BandpassFilter(signal, 1000, [150 250]);
% 
% % get ripple duration in sec
% ripple_dur = ripples(:,3) - ripples(:,1);
% 
% % pull out unfiltered and filtered ripples from lfp
% for r = 1:size(ripples,1)
%     idx = lfp_ts >= ripples(r,1) & lfp_ts <= ripples(r,3);
%     unfiltered_ripple{r} = signal(idx);
%     filtered_ripple{r} = signal_filtered(idx);
%     time{r} = lfp_ts(idx);
% end

% pack into struct
for ch = 1:length(all_ripples)
    ripple_info.ripple_start{ch} = all_ripples{ch}(:,1);
    ripple_info.ripple_peak{ch} = all_ripples{ch}(:,2);
    ripple_info.ripple_end{ch} = all_ripples{ch}(:,3);
    ripple_info.peakNormalizedPower{ch} = all_ripples{ch}(:,4);
    ripple_info.frequency{ch} = all_ripples{ch}(:,5);
    ripple_info.amplitude{ch} = all_ripples{ch}(:,6);
    ripple_info.ripple_dur{ch} = ripple_dur{ch};
    ripple_info.unfiltered_ripple{ch} = unfiltered_ripple{ch};
    ripple_info.filtered_ripple{ch} = filtered_ripple{ch};
    ripple_info.time{ch} = time{ch};
end


if figs
    for ch = 1:length(all_ripples)
        figure;
        p = ceil(sqrt(size(all_ripples{ch},1)));
        colors = viridis(length(unique(all_ripples{ch}(:,4))));
        for r = 1:size(all_ripples{ch},1)
            subplot(p,p,r)
            plot(ripple_info.time{ch}{r},ripple_info.unfiltered_ripple{ch}{r},...
                'color',interp1(unique(all_ripples{ch}(:,4)),colors,all_ripples{ch}(r,4)))
            hold on
            box off
            axis off
        end
        sgtitle([num2str(r),' Ripples on Channel ',num2str(ch),', Color code: power'],'Color','w')
        darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
    end
    
    
    
    
%     figure;
%     p = ceil(sqrt(size(ripples,1)));
%     colors = viridis(length(unique(ripples(:,4))));
%     for r = 1:size(ripples,1)
%         subplot(p,p,r)
%         plot(time{r},unfiltered_ripple{r},'color',interp1(unique(ripples(:,4)),colors,ripples(r,4)))
%         hold on
%         box off
%         axis off
%     end
%     sgtitle([num2str(r),' Ripples on Channel ',num2str(ch),', Color code: power'],'Color','w')
%     darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
    
%     figure
%     for r = 1:size(ripples,1)
%         subplot(p,p,r)
%         idx = lfp_ts >= ripples(r,1)-0.5 & lfp_ts <= ripples(r,3)+0.5;
%         pspectrum(signal(idx),frequency,'spectrogram',...
%             'OverlapPercent',99,'MinTHreshold',10,'FrequencyLimits',[100, 250])
%         colorbar off
%         axis off
%         box off
%         title('')
%         colormap(viridis)
%         hold on
%         pause(.0001)
%     end
    
%     figure;
%     plot(lfp_ts,signal);
%     hold on
%     yLim = ylim;
%     for j=1:size(ripples,1)
%         plot([ripples(j,1) ripples(j,1)],yLim,'g-');
%         plot([ripples(j,2) ripples(j,2)],yLim,'k-');
%         plot([ripples(j,3) ripples(j,3)],yLim,'r-');
%     end
end
end


%FindRipples - Find hippocampal ripples (~200Hz oscillations).
% 
%  USAGE
%
%    ripples = FindRipples(samples,<options>)
%
%    samples        filtered LFP (one channel).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'show'        plot results (default = 'off')
%     'noise'       noisy ripple-band filtered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%    =========================================================================
%
%  OUTPUT
%
%    ripples        for each ripple, [start_t peak_t end_t peakNormalizedPower fqcy]
%

% Copyright (C) 2004-2006 by Michaël Zugaro, adapted from Hajime Hirase, modified by A Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function [ripples,sd] = findRipples(samples,varargin)

% Default values
frequency = 1250;
show = 'off';

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help FindRipples'' for details).');
end

% Parameters
windowLength = floor(frequency/1250*11);
minRippleLength = 25;
maxRippleLength = 350;
minInterRippleInterval = 30;
% Ripple envoloppe must exceed lowThresholdFactor*stdev
lowThresholdFactor = 3;
% Ripple peak must exceed highThresholdFactor*stdev
highThresholdFactor = 7;
%  % Check parameter sizes
%  if size(samples,2) ~= 2,
%  	error('Parameter ''samples'' is not a Nx2 matrix (type ''help FindRipples'' for details).');
%  end
sd = 0;

% Parse parameter list
for i = 1:2:length(varargin),
	if ~isa(varargin{i},'char'),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help FindRipples'' for details).']);
	end
	switch(lower(varargin{i})),
        case 'thresholds',
			thresholds = varargin{i+1};
% 			if ~isivector(thresholds,'#2','>0'),
% 				error('Incorrect value for property ''thresholds'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
% 			end
			lowThresholdFactor = thresholds(1);
			highThresholdFactor = thresholds(2);
		case 'frequency',
			frequency = varargin{i+1};
			if ~isa(frequency,'numeric') | length(frequency) ~= 1 | frequency <= 0,
				error('Incorrect value for property ''frequency'' (type ''help FindRipples'' for details).');
			end
		case 'show',
			show = varargin{i+1};
%  			if ~IsStringInList(show,'on','off'),
%  				error('Incorrect value for property ''show'' (type ''help FindRipples'' for details).');
%  			end

        case 'sd',
			sd = varargin{i+1};
% 			if size(noise,1) ~= size(samples,1)
% 				error('Incorrect value for property ''noise'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
% 			end
		otherwise,
% 			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help FindRipples'' for details).']);
	end
end

ripples = [];

% Square and normalize signal
signal = samples;

squaredSignal = signal.^2;
window = ones(windowLength,1)/windowLength;

normalizedSquaredSignal = Filter0(window,sum(squaredSignal,2));

if ~sd
    sd = std(normalizedSquaredSignal);
end
normalizedSquaredSignal = unity(normalizedSquaredSignal,sd);


% Detect ripple periods by thresholding normalized squared signal
thresholded = normalizedSquaredSignal > lowThresholdFactor;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
if length(stop) == length(start)-1,
	% Exclude last ripple if it is incomplete
	start = start(1:end-1);
end
if length(stop)-1 == length(start)
	% Exclude first ripple if it is incomplete
    stop = stop(2:end);
end
firstPass = [start,stop];
	disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);

%added by A Peyrache - exclude ripples whos length<minRippleLength

minIx = stop - start < minRippleLength/1000*frequency;
stop(minIx) = [];
start(minIx) = [];

minIx = stop - start > maxRippleLength/1000*frequency;
stop(minIx) = [];
start(minIx) = [];

firstPass = [start,stop];

if isempty(firstPass),
	disp('Detection by thresholding failed');
	return
else
	disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

% Merge ripples if inter-ripple period is too short
minInterRippleSamples = minInterRippleInterval/1000*frequency;
secondPass = [];
ripple = firstPass(1,:);
for i = 2:size(firstPass,1)
	if firstPass(i,1) - ripple(2) < minInterRippleSamples,
		% Merge
		ripple = [ripple(1) firstPass(i,2)];
	else
		secondPass = [secondPass ; ripple];
		ripple = firstPass(i,:);
	end    
end
secondPass = [secondPass ; ripple];
if isempty(secondPass),
	disp('Ripple merge failed');
    ripples = [];
	return
else
	disp(['After ripple merge: ' num2str(length(secondPass)) ' events.']);
end

% Discard ripples with a peak power < highThresholdFactor
thirdPass = [];
peakNormalizedPower = [];
for i = 1:size(secondPass,1)
	[maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
	if maxValue > highThresholdFactor,
		thirdPass = [thirdPass ; secondPass(i,:)];
		peakNormalizedPower = [peakNormalizedPower ; maxValue];
	end
end
if isempty(thirdPass),
	disp('Peak thresholding failed.');
	return
else
	disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end

% Detect negative peak position for each ripple
peakPosition = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1),
	[minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
	peakPosition(i) = minIndex + thirdPass(i,1) - 2;
end

% Average Fqcy (by means of average inter-peak intervals) Added  by A
% Peyrache

fqcy = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1),
	peakIx = LocalMinima(signal(thirdPass(i,1):thirdPass(i,2)),4,0);
    if ~isempty(peakIx)
        fqcy(i) = frequency/median(diff(peakIx));
    else
        warning('Weird...')
%         keyboard
    end
end

% answer of this function...
ripples = [thirdPass(:,1)/frequency peakPosition/frequency thirdPass(:,2)/frequency peakNormalizedPower fqcy];

% Optionally, plot results
if strcmp(show,'on'),
	figure;
	time = (0:length(signal)-1)/frequency;
	subplot(3,1,1);hold on;
	plot(time,signal);
	subplot(3,1,2);hold on;
	plot(time,squaredSignal);
	subplot(3,1,3);hold on;
	plot(time,normalizedSquaredSignal);
	for i = 1:3,
		subplot(3,1,i);
		yLim = ylim;
		for j=1:size(ripples,1),
			plot([ripples(j,1) ripples(j,1)],yLim,'g-');
			plot([ripples(j,2) ripples(j,2)],yLim,'k-');
			plot([ripples(j,3) ripples(j,3)],yLim,'r-');
			if i == 3,
				plot([ripples(j,1) ripples(j,3)],[ripples(j,4) ripples(j,4)],'k-');
			end
		end
	end
end
end

function y = Filter0(b, x)

if size(x,1) == 1
	x = x(:);
end

if mod(length(b),2)~=1
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];
end

function U = unity(A,stdA)

meanA = mean(A);
U = (A - repmat(meanA,size(A,1),1))./repmat(stdA,size(A,1),1);

end


% mins = LocalMinima(x, NotCloserThan, LessThan)
%
% finds positions of all local minima in input array
% in the case that it goes down, stays level, then goes up,
% the function finds the earliest point of equality
%
% second optional argument gives minimum distance between
% two minima - if 2 are closer than this, it choses the lower between them.
%
% 3rd optional argument lets you only have minima less than a certain number
% use this option if you are computing minima of a long array - it'll take way
% less time and memory.
%
% endpoints will not be counted as minima.

% This program is the curse of my life.  Why can't things be simple?

function Mins = LocalMinima(x, NotCloserThan, LessThan)

if min(size(x))>1
    error('x should be a vector');
end
x = x(:);
nPoints = length(x);


% only look at those below LessThan
if nargin<3
    BelowThresh = (1:nPoints)';
else
    BelowThresh = find(x<LessThan);
end
xBelow = x(BelowThresh);
GapToLeft = find(diff([0; BelowThresh])>1);
GapToRight = find(diff([BelowThresh; nPoints+1])>1);

% compute left and right signs, bearing in mind that some points are missing
sDiff = sign(diff(xBelow));
LeftSign = [1; sDiff];
LeftSign(GapToLeft) = -1;
RightSign = [sDiff; -1];
RightSign(GapToRight) = 1;

% OK, now any zero right signs need to be replaced with next non-zero ...
Zeros = find(RightSign==0);
for i=fliplr(Zeros(:)')
    RightSign(i) = RightSign(i+1);
end
    
% now we can find local minima
Mins = BelowThresh(find(LeftSign<0 & RightSign>0));

% now remove any that are too close

if nargin>=2
    while 1
        TooClose = find(diff(Mins)<NotCloserThan);
        if isempty(TooClose)
            break;
        end
%        Vals = x(Mins(TooClose:TooClose+1));
        Vals = [x(Mins(TooClose)) , x(Mins(TooClose+1))];
        [dummy Offset] = max(Vals,[],2);
        Delete = TooClose + Offset -1;
        Mins(unique(Delete)) = [];
    end
end

return

if nargin<3
    LessThan = inf;
end



nPoints = length(x);
nMins = 0;
ArraySize = floor(nPoints/10);
Mins = zeros(ArraySize,1);
PrevSign = 1;
for i=1:length(x)-1
    NextSign = sign(x(i+1)-x(i));
    
    % do we have a minimum?
    if (PrevSign<0 & NextSign>0 & x(i)<LessThan)
        nMins = nMins+1;
        Mins(nMins) = i;
    end

    % reset PrevSign, if we are not in equality situation
    if NextSign
        PrevSign=NextSign;
    end
end

% use only those we have
if nMins<ArraySize
    Mins(nMins+1:ArraySize) = [];
end
    
% look for duplicates    

if nargin>=2
    while 1
        TooClose = find(diff(Mins)<NotCloserThan);
        if isempty(TooClose)
            break;
        end
        Vals = x(Mins(TooClose:TooClose+1));
        [dummy Offset] = max(Vals,[],2);
        Delete = TooClose + Offset -1;
        Mins(unique(Delete)) = [];
    end
end

return


nPoints = length(x);

% use a trick to save memory - with findstr
s = int8([1 sign(diff(x))]);


Zeros = find(s==0);
NonZeros = uint32(find(s~=0));

% wipe zeros ... 
s(Zeros) = [];

%mins = find(s(1:nPoints-1)==-1 & s(2:end)==1);
mins = double(NonZeros(findstr(s, [-1 1])));

if nargin>=3
    mins = mins(find(x(mins)<LessThan));
end

if nargin>=2
    while 1
        TooClose = find(diff(mins)<NotCloserThan);
        if isempty(TooClose)
            break;
        end
        Vals = x(mins(TooClose:TooClose+1));
        [dummy Offset] = max(Vals,[],2);
        Delete = TooClose + Offset -1;
        mins(unique(Delete)) = [];
    end
end
end