%PlotPopEPSP - Plot population EPSPs (waveforms, slope, amplitude) across time.
%
%  Slope is determined as the maximum slope during the time window around stimulation.
%  Amplitude is determined as the maximum amplitude during the same time window.
%
%  USAGE
%
%    stats = PlotPopEPSP(stims,lfp,<options>)
%
%    stims          stimulation timestamps
%    lfp            local field potential samples
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'direction'   direction of the population EPSP, either 'up' (recording
%                   from str. pyr., default) or 'down' (recording from str.
%                   radiatum)
%     'durations'   durations before and after synchronizing events for each
%                   trial (in s) (default = [-0.01 0.03])
%     'focus'       window (in s) where max amplitude and slope should be
%                   sought (default = [0.003 0.015])
%     'window'      window length (in s or in counts) for running averages
%                   (default = 100)
%     'mode'        compute slope and amplitude running averages using a
%                   fixed time window ('time'), or a window with a fixed
%                   number of stimulations ('count', default)
%     'parent'      parent figure or uipanel handle (default = gcf)
%     'plot'        either 'scatter' (default) or 'average'
%     'debug'       plot debugging information (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    stats.amplitude.beta      coefficients of linear regression for amplitude
%                              vs time (or count)
%    stats.amplitude.p         p-values for the above coefficients
%    stats.amplitude.change    percent increase or decrease per second (or count)
%    stats.slope.beta          coefficients for slope vs time (or count)
%    stats.slope.p             p-values for the above coefficients
%    stats.slope.change        percent increase or decrease per second (or count)

% Copyright (C) 2010-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function stats = PlotPopEPSP(stims,lfp,varargin)

% Default values
parent = [];
durations = [-0.01 0.03];
mode = 'count';
window = 100;
focus = [0.003 0.015];
debug = 'off';
what = 'scatter';
direction = 'up';
stats = [];

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
end

% Check parameter sizes
if ~isdvector(stims),
	error('Parameter ''stims'' is not a vector (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
end
if ~isdmatrix(lfp),
	error('Parameter ''lfp'' is not a matrix (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'direction',
			direction = lower(varargin{i+1});
			if ~isstring(direction,'up','down'),
				error('Incorrect value for property ''direction'' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
			end
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','<'),
				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
			end
		case 'focus',
			focus = varargin{i+1};
			if ~isdvector(focus,'#2','<'),
				error('Incorrect value for property ''focus'' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
			end
		case 'mode',
			mode = lower(varargin{i+1});
			if ~isstring(mode,'time','count'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
			end
		case 'parent',
			parent = varargin{i+1};
			if ~ishandle(parent),
				error('Incorrect value for property ''parent'' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
			end
		case 'plot',
			what = lower(varargin{i+1});
			if ~isstring(what,'scatter','average'),
				error('Incorrect value for property ''what'' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
			end
		case 'debug',
			debug = lower(varargin{i+1});
			if ~isstring(debug,'on','off'),
				error('Incorrect value for property ''debug'' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotPopEPSP">PlotPopEPSP</a>'' for details).']);
	end
end

up = strcmp(direction,'up');
if isempty(parent), parent = gcf; end

% Some analyses depend on the method (fixed time vs fixed count), but others do not
% To simplify the code, we will use two variables:
%  - stims0 uses time, whatever the method
%  - stims  uses either time or count, depending on the method
stims0 = stims;
if strcmp(mode,'count'),
	stims = (1:length(stims))';
end

% Constants
nStims = length(stims);
height = max(lfp(:,2))-min(lfp(:,2));
nSlices = 10;
start = stims(1);
stop = stims(end);
window = (stop-start)/nSlices;
slices = (0:nSlices-1)*window+start;

mainFigure = gcf;
if strcmp(debug,'on'),
	debugFigure = figure;
	figure(mainFigure);
end
panel1 = Subpanel(2,3,[1 4],'parent',parent);
panel2 = Subpanel(2,3,[2 3],'parent',parent);
panel3 = Subpanel(2,3,[5 6],'parent',parent);

% 1) Process pop EPSPs in each time (or count) slice

for i = 1:nSlices,

	% Time slice
	ok = InIntervals(stims,[0 window]+start+(i-1)*window);
	if ~any(ok), continue; end

	% Resynchronize pop EPSPs, compute mean and conf intervals, and distributions
	[sync,index] = Sync(lfp,stims0(ok),'durations',durations);
	[average,unused,conf,time] = SyncHist(sync,index,'mode','mean','durations',durations);
	[d,x] = SyncHist(sync,index,'mode','dist','durations',durations,'smooth',[1 1]);
	% Store distributions in a larger matrix (this will contain distributions for all time slices)
	n = size(d,1);
	dist((i-1)*n+(1:n),:) = d;

	% Plot individual traces + mean and conf intervals
	% (distributions are not plotted now, they will be plotted all at once at the end)
	subplot(1,2,1,'Parent',panel1);
	shift = (i-1)*height*0.75;
	PlotSync([sync(:,1) sync(:,2)+shift],index,'durations',durations);
	PlotMean(time,average+shift,conf(:,1)+shift,conf(:,2)+shift,'-','k');

	if i == 1, m = min(sync(:,2)); elseif i == nSlices, M = max(sync(:,2)+shift); end

end

% Adjust axes for individual traces + mean and conf intervals
subplot(1,2,1,'Parent',panel1);
ylim([m M]);
set(gca,'ytick',(0:nSlices-1)*height*0.75,'yticklabel',slices);
xlabel('Time (s)');
ylabel('Population EPSP');

% Plot distributions and adjust axes
subplot(1,2,2,'Parent',panel1);
PlotColorMap(dist,1,'x',x);
ylim([0 nSlices*n]);
set(gca,'ytick',[]);
xlabel('Time (s)');


% 2) Compute time-varying slope and amplitude

% Resynchronize all pop EPSPs (use only 15ms after stims)
[sync,index] = Sync(lfp,stims0,'durations',[-0.005 0.015]);

% Compute slope and amplitude for each stimulation
figure(debugFigure);
nStims = max(index);
frequency = 1/median(diff(lfp(:,1)));
smooth = frequency/1250;
nStimsPlottedPerSlice = 5;
nStimsPerSlice = nStims/nSlices;
nSkip = floor(nStims/(nSlices*(nStimsPlottedPerSlice-1))); % plot one stim every n
colors = Bright(nStimsPlottedPerSlice);
for i = 1:nStims,
	% Select the appropriate data, smooth and differentiate
	s = sync(index==i,:);
	s(:,2) = Smooth(s(:,2),smooth);
	ds = Diff(s);
	% To compute the amplitude, we need to vertically align all traces to the baseline level,
	% i.e. the level at approx.  t=0
	atStim = s(:,1)>-0.001&s(:,1)<0.001;
	m = mean(s(atStim,2));
	s(:,2) = s(:,2) - m;
	% Compute amplitude as the first local maximum (or minimum) in the appropriate time interval
	in = InIntervals(s(:,1),focus);
	if up,
		extrema = find(IsExtremum(s)&in);
	else
		extrema = find(IsExtremum(s,'mode','minima')&in);
	end
	if isempty(extrema),
		amplitude(i) = nan;
		tAmplitude(i) = nan;
		slope(i) = nan;
		tSlope(i) = nan;
		continue;
	end
	extremum = extrema(1);
	tAmplitude(i) = s(extremum,1);
	amplitude(i) = s(extremum,2);
	% Compute slope as the maximum slope in the appropriate time interval and before the max amplitude
	% (set all values outside this time interval to -inf to make sure the max will be in the interval)
	dS = ds;
	dS(~in,2) = nan;
	afterAmplitude = dS(:,1)>tAmplitude(i);
	dS(afterAmplitude,2) = nan;
	if up,
		[slope(i),j] = nanmax(dS(:,2));
	else
		[slope(i),j] = nanmin(dS(:,2));
	end
	tSlope(i) = dS(j,1);
	% Debug: plot slopes and amplitudes
	slice = floor(i/nStimsPerSlice);
	stimInSlice = floor(i-1-slice*nStimsPerSlice);
	if strcmp(debug,'on') && rem(stimInSlice,nSkip) == 0,
		k = ceil(nSlices*i/nStims);
		l = floor(stimInSlice/nSkip)+1;
		SquareSubplot(nSlices,k);
		xlabel(num2str(slices(k)));
		hold on;
		PlotXY(s,'color',colors(l,:),'linewidth',3);
		PlotSlope(tAmplitude(i),amplitude(i),0,0.005,'r');
		PlotSlope(tSlope(i),s(j,2),slope(i),0.005,'k');
	end
end
if strcmp(debug,'on'),
	% Make sure all subplots use the same scale
	sub = get(debugFigure,'children');
	if ~isempty(sub),
		for i = 1:length(sub),
			lims(i,:) = ylim(sub(i));
		end
		m = min(lims(:,1));
		M = max(lims(:,2));
		for i = 1:length(sub),
			ylim(sub(i),[m M]);
		end
	end
end

% For linear regressions and running averages, 'rectify' slope and amplitude
if ~up,
	amplitude = -amplitude;
	slope = -slope;
end

% Compute and plot scatterplots or running averages
figure(mainFigure);
if strcmp(mode,'time'),
	xLabel = 'Time (s)';
else
	xLabel = 'Count';
end
if strcmp(what,'scatter'),
	% Slope (scatterplot)
	good = ~isnan(slope);
	st = regstats(slope(good),stims(good),'linear','tstat');
	stats.slope.beta = st.tstat.beta;
	stats.slope.change = 100*stats.slope.beta(2)/stats.slope.beta(1);
	stats.slope.p = st.tstat.pval;
	if stats.slope.p < 0.05, color = 'r'; else color = 'k'; end
	subplot(2,3,1,'Parent',panel2);
	hold on;
	plot(stims,slope,'.');
	x0 = (stims(1,1)+stims(end,1))/2;
	dx = stims(end,1)-stims(1,1);
	y0 = stats.slope.beta(2)*x0+stats.slope.beta(1);
	PlotSlope(x0,y0,stats.slope.beta(2),dx,color,'linewidth',2);
	xlim([stims(1,1) stims(end,1)]);
	ylabel('Slope (a.u.)');
	set(gca,'xtick',[]);
	% Slope time (scatterplot)
	good = ~isnan(tSlope);
	st = regstats(tSlope(good),stims(good),'linear','tstat');
	betaTSlope = st.tstat.beta;
	pTSlope = st.tstat.pval;
	if pTSlope < 0.05, color = 'r'; else color = 'k'; end
	subplot(2,3,2,'Parent',panel2);
	hold on;
	plot(stims,tSlope,'.');
	x0 = (stims(1,1)+stims(end,1))/2;
	dx = stims(end,1)-stims(1,1);
	y0 = betaTSlope(2)*x0+betaTSlope(1);
	PlotSlope(x0,y0,betaTSlope(2),dx,color,'linewidth',2);
	xlim([stims(1,1) stims(end,1)]);
	ylabel('Slope time (ms)');
	set(gca,'xtick',[]);
else
	% Slope (running average)
	[ts,ms,es] = RunningAverage(stims,slope,'window',window,'overlap',0.5*window);
	subplot(2,3,1,'Parent',panel2);
	PlotMean(ts,ms,es(:,1),es(:,2),':','k');
	xlim([stims(1,1) stims(end,1)]);
	PlotHVLines(slices,'v','color',[0.8 0.8 0.8]);
	ylabel('Slope (a.u.)');
	set(gca,'xtick',[]);
	% Slope time (running average)
	[ts,ms,es] = RunningAverage(stims,tSlope*1000,'window',window,'overlap',0.5*window);
	subplot(2,3,2,'Parent',panel2);
	PlotMean(ts,ms,es(:,1),es(:,2),':','k');
	xlim([stims(1,1) stims(end,1)]);
	PlotHVLines(slices,'v','color',[0.8 0.8 0.8]);
	ylabel('Slope time (ms)');
	set(gca,'xtick',[]);
end
if strcmp(what,'scatter'),
	% Amplitude (scatter)
	st = regstats(amplitude(good),stims(good),'linear','tstat');
	stats.amplitude.beta = st.tstat.beta;
	stats.amplitude.change = 100*stats.amplitude.beta(2)/stats.amplitude.beta(1);
	stats.amplitude.p = st.tstat.pval;
	if stats.amplitude.p < 0.05, color = 'r'; else color = 'k'; end
	subplot(2,3,4,'Parent',panel2);
	hold on;
	plot(stims,amplitude,'.');
	y0 = stats.amplitude.beta(2)*x0+stats.amplitude.beta(1);
	PlotSlope(x0,y0,stats.amplitude.beta(2),dx,color,'linewidth',2);
	xlim([stims(1,1) stims(end,1)]);
	xlabel(xLabel);
	ylabel('Amplitude (a.u.)');
	% Amplitude time (scatter)
	st = regstats(tAmplitude(good),stims(good),'linear','tstat');
	betaTAmplitude = st.tstat.beta;
	pTAmplitude = st.tstat.pval;
	if pTAmplitude < 0.05, color = 'r'; else color = 'k'; end
	subplot(2,3,5,'Parent',panel2);
	hold on;
	plot(stims,tAmplitude,'.');
	y0 = betaTAmplitude(2)*x0+betaTAmplitude(1);
	PlotSlope(x0,y0,betaTAmplitude(2),dx,color,'linewidth',2);
	xlim([stims(1,1) stims(end,1)]);
	xlabel(xLabel);
	ylabel('Amplitude time (ms)');
else
	% Amplitude (running average)
	[ta,ma,ea] = RunningAverage(stims,amplitude,'window',window,'overlap',0.5*window);
	subplot(2,3,4,'Parent',panel2);
	PlotMean(ta,ma,ea(:,1),ea(:,2),':','k');
	xlim([stims(1,1) stims(end,1)]);
	PlotHVLines(slices,'v','color',[0.8 0.8 0.8]);
	xlabel(xLabel);
	ylabel('Amplitude (a.u.)');
	% Amplitude time (running average)
	[ta,ma,ea] = RunningAverage(stims,tAmplitude*1000,'window',window,'overlap',0.5*window);
	subplot(2,3,5,'Parent',panel2);
	PlotMean(ta,ma,ea(:,1),ea(:,2),':','k');
	xlim([stims(1,1) stims(end,1)]);
	PlotHVLines(slices,'v','color',[0.8 0.8 0.8]);
	ylabel('Amplitude time (ms)');
	xlabel(xLabel);
end

% 3) CV

binSize = 0.1;
smooth = 10;
for i = 1:20,
	cv(i) = CV(stims0,'order',i,'binSize',binSize,'smooth',smooth,'method','fixed');
end
cv2 = CV(stims0,'measure','cv2');
subplot(2,3,3,'Parent',panel2);
% ISIs
ds = diff(stims0(:,1));
x = 0:0.1:10;
h = hist(ds,x);
h = h/sum(h);
bar(x,h);
xlabel('ISI');
xlim([0 10]);
ylabel('Count');
% CV and CV2
subplot(2,3,6,'Parent',panel2);
plot(cv,'+-');
hold on;
plot(cv2,'r+');
xlabel('nth ISI');
ylabel('CV & CV2');

% 4) Stims autocorrelograms

[ccg,x,y] = ShortTimeCCG(stims0,'min',100,'mode','norm','smooth',[2 0]);
PlotShortTimeCCG(ccg,'x',x,'y',y,'Parent',panel3);
clim([0 0.01]);
