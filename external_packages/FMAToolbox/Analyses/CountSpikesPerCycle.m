%CountSpikesPerCycle - Count number of spikes per LFP cycle.
%
% Count the number of spikes per cycle in the ongoing oscillatory LFP
% (e.g. during theta).
%
%  USAGE
%
%    [count,used] = CountSpikesPerCycle(spikes,phases)
%
%    spikes         list of spike timestamps
%    phases         instantaneous phases in radians (see <a href="matlab:help Phase">Phase</a>)
%
%  OUTPUT
%
%    count          spike counts in successive cycles
%    used           logical vector indicating which spikes were used
%                   (and which were discarded because they belonged
%                   to incomplete cycles)
%
%  SEE
%
%    See also Phase.

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [count,used] = CountSpikesPerCycle(spikes,phases)

% Check number of parameters
if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help CountSpikesPerCycle">CountSpikesPerCycle</a>'' for details).');
end

% Check parameter sizes
if ~isdvector(spikes),
	error('Parameter ''spikes'' is not a vector (type ''help <a href="matlab:help CountSpikesPerCycle">CountSpikesPerCycle</a>'' for details).');
end
if ~isdmatrix(phases) | size(phases,2) ~= 2,
	error('Parameter ''phases'' is not an Nx2 matrix (type ''help <a href="matlab:help CountSpikesPerCycle">CountSpikesPerCycle</a>'' for details).');
end
isradians(phases(:,2));

% Find theta peaks
p = phases(:,2);
p(p>pi) = p(p>pi) - 2*pi;
[up,unused] = ZeroCrossings([phases(:,1) p]);
cycles = find(up);

% Intervals between successive theta peaks
[unused,interval] = InIntervals(spikes,[phases(cycles(1:end-1),1) phases(cycles(2:end),1)]);

% Count
used = interval ~= 0;
count = Accumulate(interval(used));
count = count(interval(used));

