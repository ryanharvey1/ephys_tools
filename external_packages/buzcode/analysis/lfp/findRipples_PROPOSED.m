function [ripples,sd,bad] = bz_FindRipples(varargin)
%FindRipples - Find hippocampal ripples (100~200Hz oscillations).
%
%  USAGE
%
%    [ripples,stdev,noise] = bz_FindRipples(filtered,timestamps,<options>)
%
%    OR
%
%   [ripples,stdev,noise] = bz_FindRipples(baspath,channel,<options>)
%
%    Ripples are detected using the normalized squared signal (NSS) by
%    thresholding the baseline, merging neighboring events, thresholding
%    the peaks, and discarding events with excessive duration.
%    Thresholds are computed as multiples of the standard deviation of
%    the NSS. Alternatively, one can use explicit values, typically obtained
%    from a previous call.
%
%    filtered       ripple-band filtered LFP (one channel).
%	 timestamps	    timestamps to match filtered variable
%    <options>      optional list of property-value pairs (see table below)
%
% %%%%  OR BELOW    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%
%    basepath       path to a single session to run findRipples on
%    channel      	Ripple channel to use for detection
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'thresholds'  thresholds for ripple beginning/end and peak, in multiples
%                   of the stdev (default = [2 5])
%     'durations'   min inter-ripple interval and max ripple duration, in ms
%                   (default = [30 100])
%     'baseline'    interval used to compute normalization (default = all)
%     'restrict'    same as 'baseline' (for backwards compatibility)
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'stdev'       reuse previously computed stdev
%     'show'        plot results (default = 'off')
%     'noise'       noisy ripple-band filtered channel used to exclude ripple-
%                   like noise (events also present on this channel are
%                   discarded)
%    =========================================================================
%
%  OUTPUT
%
%    ripples        for each ripple, [start_t peak_t end_t peakNormalizedPower]
%    stdev          standard deviation of the NSS (can be reused subsequently)
%    noise          ripple-like activity recorded simultaneously on the noise
%                   channel (for debugging info)
%
%  SEE
%
%    See also FilterLFP, RippleStats, SaveRippleEvents, PlotRippleStats.

% Copyright (C) 2004-2011 by Michaël Zugaro, initial algorithm by Hajime Hirase
% edited by David Tingley, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


warning('this function is under development and may not work... yet')

% Default values
p = inputParser;
addParameter(p,'frequency',1250,@isnumeric)
addParameter(p,'show','off',@isstr)
addParameter(p,'thresholds',[2 5],@isivector)
addParameter(p,'durations',[30 100],@isnumberic)
addParameter(p,'restrict',[],@isnumeric)
addParameter(p,'stdev',[],@isnumeric)
addParameter(p,'noise',[],@isdmatrix)

if isstr(varargin{1})  % if first arg is basepath
    addRequired(p,'basepath',@isstr)
    addRequired(p,'channel',@isnumeric)    
    parse(p,varargin{:})
    
    xmlfile = dir([p.Results.basepath '/*xml']);
    SetCurrentSession([p.Results.basepath '/' xmlfile.name]);
    lfp = GetLFP(p.Results.channel);
    
    filtered = FilterLFP(lfp,'passband',[120 220]);
elseif isnumeric(varargin{1}) % if first arg is filtered LFP
    addRequired(p,'filtered',@isnumeric)
    addRequired(p,'timestamps',@isnumeric)
    parse(p,varargin{:})
    filtered = p.Results.filtered;
    timestamps = p.Results.timestamps;
end

% assign parameters (either defaults or given)
frequency = p.Results.frequency;
show = p.Results.show;
restrict = p.Results.restrict;
sd = p.Results.stdev;
noise = p.Results.noise;
lowThresholdFactor = p.Results.thresholds(1);
highThresholdFactor = p.Results.thresholds(2);
minInterRippleInterval = p.Results.durations(1);
maxRippleDuration = p.Results.durations(2);

% Check parameter sizes
if ~isdmatrix(filtered) | size(filtered,2) ~= 2,
	error('Parameter ''filtered'' is not a Nx2 matrix (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
end

% Parameters
windowLength = frequency/1250*11;

% Square and normalize signal
signal = filtered(:,2);
squaredSignal = signal.^2;
window = ones(windowLength,1)/windowLength;
keep = [];
if ~isempty(restrict),
	keep = filtered(:,1)>=restrict(1)&filtered(:,1)<=restrict(2);
end

[normalizedSquaredSignal,sd] = unity(Filter0(window,sum(squaredSignal,2)),sd,keep);

% Detect ripple periods by thresholding normalized squared signal
thresholded = normalizedSquaredSignal > lowThresholdFactor;
start = find(diff(thresholded)>0);
stop = find(diff(thresholded)<0);
% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1,
	start = start(1:end-1);
end
% Exclude first ripple if it is incomplete
if length(stop)-1 == length(start),
    stop = stop(2:end);
end
% Correct special case when both first and last ripples are incomplete
if start(1) > stop(1),
	stop(1) = [];
	start(end) = [];
end
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
	peakPosition(i) = minIndex + thirdPass(i,1) - 1;
end

% Discard ripples that are way too long
time = filtered(:,1);
ripples = [time(thirdPass(:,1)) time(peakPosition) time(thirdPass(:,2)) peakNormalizedPower];
duration = ripples(:,3)-ripples(:,1);
ripples(duration>maxRippleDuration/1000,:) = [];
disp(['After duration test: ' num2str(size(ripples,1)) ' events.']);

% If a noisy channel was provided, find ripple-like events and exclude them
bad = [];
if ~isempty(noise),
	% Square and pseudo-normalize (divide by signal stdev) noise
	squaredNoise = noise(:,2).^2;
	window = ones(windowLength,1)/windowLength;
	normalizedSquaredNoise = unity(Filter0(window,sum(squaredNoise,2)),sd,[]);
	excluded = logical(zeros(size(ripples,1),1));
	% Exclude ripples when concomittent noise crosses high detection threshold
	previous = 1;
	for i = 1:size(ripples,1),
		j = FindInInterval(noise,[ripples(i,1),ripples(i,3)],previous);
		previous = j(2);
		if any(normalizedSquaredNoise(j(1):j(2))>highThresholdFactor),
			excluded(i) = 1;
		end
	end
	bad = ripples(excluded,:);
	ripples = ripples(~excluded,:);
	disp(['After noise removal: ' num2str(size(ripples,1)) ' events.']);
end

% Optionally, plot results
if strcmp(show,'on'),
	figure;
	if ~isempty(noise),
		MultiPlotXY([time signal],[time squaredSignal],[time normalizedSquaredSignal],[time noise(:,2)],[time squaredNoise],[time normalizedSquaredNoise]);
		nPlots = 6;
		subplot(nPlots,1,3);
 		ylim([0 highThresholdFactor*1.1]);
		subplot(nPlots,1,6);
  		ylim([0 highThresholdFactor*1.1]);
	else
		MultiPlotXY([time signal],[time squaredSignal],[time normalizedSquaredSignal]);
%  		MultiPlotXY(time,signal,time,squaredSignal,time,normalizedSquaredSignal);
		nPlots = 3;
		subplot(nPlots,1,3);
  		ylim([0 highThresholdFactor*1.1]);
	end
	for i = 1:nPlots,
		subplot(nPlots,1,i);
		hold on;
  		yLim = ylim;
		for j=1:size(ripples,1),
			plot([ripples(j,1) ripples(j,1)],yLim,'g-');
			plot([ripples(j,2) ripples(j,2)],yLim,'k-');
			plot([ripples(j,3) ripples(j,3)],yLim,'r-');
			if i == 3,
				plot([ripples(j,1) ripples(j,3)],[ripples(j,4) ripples(j,4)],'k-');
			end
		end
		for j=1:size(bad,1),
			plot([bad(j,1) bad(j,1)],yLim,'k-');
			plot([bad(j,2) bad(j,2)],yLim,'k-');
			plot([bad(j,3) bad(j,3)],yLim,'k-');
			if i == 3,
				plot([bad(j,1) bad(j,3)],[bad(j,4) bad(j,4)],'k-');
			end
		end
		if mod(i,3) == 0,
			plot(xlim,[lowThresholdFactor lowThresholdFactor],'k','linestyle','--');
			plot(xlim,[highThresholdFactor highThresholdFactor],'k-');
		end
	end
end

function y = Filter0(b,x)

if size(x,1) == 1
	x = x(:);
end

if mod(length(b),2)~=1
	error('filter order should be odd');
end

shift = (length(b)-1)/2;

[y0 z] = filter(b,1,x);

y = [y0(shift+1:end,:) ; z(1:shift,:)];

function [U,stdA] = unity(A,sd,restrict)

if ~isempty(restrict),
	meanA = mean(A(restrict));
	stdA = std(A(restrict));
else
	meanA = mean(A);
	stdA = std(A);
end
if ~isempty(sd),
	stdA = sd;
end

U = (A - meanA)/stdA;

