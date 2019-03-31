%PhasePrecession - Compute spike phase.
%
% Compute spike phase precession using the methods of O'Keefe and Recce (1993),
% i.e. spike phase vs position, and Harris et al. (2002), i.e. spike phase vs
% spike rate.
%
%  USAGE
%
%    [data,stats] = PhasePrecession(positions,spikes,phases,<options>)
%
%    positions      position samples (linearized, normalized to [0,1])
%    spikes         spike timestamps
%    phases         phase samples (see <a href="matlab:help Phase">Phase</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'maxGap'      time gaps between successive position samples exceeding
%                   this threshold (e.g. undetects) will not be interpolated
%                   (default = 100 ms)
%    =========================================================================
%
%  OUTPUT
%
%    data.position  spike phase vs position:
%     .x             x coordinates
%     .t             spike times
%     .phase         spike phases (in radians)
%
%    data.rate      spike phase vs spike rate:
%     .r             spike rates
%     .t             spike times
%     .phase         spike phases (in radians)
%
%    Additional statistics computed using phase vs position data:
%
%    stats.center   center x for early, middle and late subfields
%    stats.mean     mean phase for early, middle and late subfields
%    stats.var      phase variance for early, middle and late subfields
%    stats.std      phase standard deviation for early, middle and late subfields
%    stats.conf     95% confidence intervals
%    stats.all      cell array containing all phases for each subfield
%                   (useful for population analyses)
%
%  SEE
%
%    See also Phase, PlotPhasePrecession.

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [data,stats] = PhasePrecession(positions,spikes,phases,varargin)

maxGap = 0.1;

% Check number of parameters
if nargin < 3 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
end

% Check parameter sizes
if ~isempty(positions) && size(positions,2) < 2,
	error('Parameter ''positions'' should have at least 2 columns (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
end
if size(spikes,2) ~= 1,
	error('Parameter ''spikes'' is not a vector (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
end
if size(phases,2) ~= 2,
	error('Parameter ''phases'' is not a Nx2 matrix (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
end
% isradians(phases(:,2));

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).']);
	end
	switch(lower(varargin{i})),

		case 'maxgap',
			maxGap = varargin{i+1};
			if ~isdscalar(maxGap,'>0'),
				error('Incorrect value for property ''maxGap'' (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PhasePrecession">PhasePrecession</a>'' for details).']);

	end
end

% Default values
data.position.x = [];
data.position.t = [];
data.position.phase = [];
data.rate.r = [];
data.rate.t = [];
data.rate.phase = [];
stats.center = [0 0 0];
stats.mean = [0 0 0];
stats.var = [0 0 0];
stats.std = [0 0 0];
stats.conf = [0 0 0;0 0 0];

% Compute spike phases
if isempty(spikes), return; end
spikePhases = Interpolate(phases,spikes,'trim','off','type','circular');
if isempty(spikePhases), return; end

% Interpolate positions at spike times
if ~isempty(positions),
	% Make sure positions are normalized
	if max(positions(:,2)) > 1 || min(positions(:,2)) < 0,
		positions(:,2) = ZeroToOne(positions(:,2));
		warning('Parameter ''positions'' should contain values in [0 1]. The data will now be transformed accordingly.');
	end
	[positions,ignored] = Interpolate(positions,spikes,'trim','off','maxGap',maxGap);
	data.position.t = positions(:,1);
	data.position.x = positions(:,2);
	data.position.phase = spikePhases(~ignored,2);
end

% Count spikes per cycle
%  p = phases(:,2);
%  p(p>pi) = p(p>pi) - 2*pi;
%  [up,unused] = ZeroCrossings([phases(:,1) p]);
%  cycles = find(up);
%  [unused,interval] = InIntervals(spikes,[phases(cycles(1:end-1),1) phases(cycles(2:end),1)]);
%  good = interval ~= 0;
%  spikesPerCycle = Accumulate(interval(good));
%  data.rate.r = spikesPerCycle(interval(good));
[data.rate.r,good] = CountSpikesPerCycle(spikes,phases);
data.rate.t = spikes(good);
data.rate.phase = spikePhases(good,2);

% Statistics: circular means and variances for early, middle and late subfields
if ~isempty(positions),
    x0 = min(data.position.x);
    x1 = max(data.position.x);
    dx = x1-x0;
    for i = 1:3,
        ok = data.position.x > x0+dx/3*(i-1) & data.position.x < x0+dx/3*i;
        if sum(ok) == 0, continue;	end
        stats.all{i} = data.position.phase(ok);
        stats.x(i) = x0+1.5*dx/3*(i-1);
        [stats.var(i),stats.std(i)] = CircularVariance(data.position.phase(ok));
        [stats.mean(i),stats.conf(:,i)] = CircularConfidenceIntervals(data.position.phase(ok));
    end
end

