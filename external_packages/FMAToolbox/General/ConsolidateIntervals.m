%ConsolidateIntervals - Consolidate intervals.
%
% Consolidate overlapping intervals, e.g. replace [10,20] [15,25] with [10,25].
%
%  USAGE
%
%    [consolidated,target] = ConsolidateIntervals(intervals,<options>)
%
%    intervals      list of intervals
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'strict'      intervals with common bounds are consolidated ('off')
%                   or kept separate ('on') (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    consolidated   consolidated intervals
%    target         for each original interval, the index of the consolidated
%                   interval to which it belongs (empty intervals yield NaN)
%
%  SEE
%
%    See also SubtractIntervals, ExcludeIntervals, InIntervals, Restrict,
%    FindInInterval, CountInIntervals, PlotIntervals.


% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [consolidated,target] = ConsolidateIntervals(intervals,varargin)

% Default values
strict = 'off';

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'strict',
			strict = lower(varargin{i+1});
			if ~isstring(strict,'on','off'),
				error('Incorrect value for property ''strict'' (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ConsolidateIntervals">ConsolidateIntervals</a>'' for details).']);
	end
end

original = intervals;
done = logical(zeros(size(intervals,1))); % to avoid retesting intervals already consolidated
if strcmp(strict,'on'),
	for i = 1:length(intervals),
		if done(i), continue; end
		% List all intervals that overlap with the current interval
		intersect = (intervals(i,1) <= intervals(:,1) & intervals(i,2) > intervals(:,1)) ...
			| (intervals(i,1) >= intervals(:,1) & intervals(i,2) <= intervals(:,2)) ...
			| (intervals(i,1) < intervals(:,2) & intervals(i,2) >= intervals(:,2));
		% Determine smallest enclosing interval
		m = min(intervals(intersect,1));
		M = max(intervals(intersect,2));
		% Consolidate
		intervals(intersect,:) = repmat([m M],sum(intersect),1);
		done(intersect) = 1;
	end
else
	for i = 1:length(intervals),
		if done(i), continue; end
		% List all intervals that overlap with the current interval
		intersect = (intervals(i,1) <= intervals(:,1) & intervals(i,2) >= intervals(:,1)) ...
			| (intervals(i,1) >= intervals(:,1) & intervals(i,2) <= intervals(:,2)) ...
			| (intervals(i,1) <= intervals(:,2) & intervals(i,2) >= intervals(:,2));
		% Determine smallest enclosing interval
		m = min(intervals(intersect,1));
		M = max(intervals(intersect,2));
		% Consolidate
		intervals(intersect,:) = repmat([m M],sum(intersect),1);
		done(intersect) = 1;
	end
end

% Sort intervals in ascending order (and store reordering information so we can reuse it later)
[intervals,order] = sortrows(intervals,1);

% Assign each consolidated interval an ID (in ascending order)
transitions = [1;find(diff(intervals(:,1))~=0)+1;length(intervals(:,1))];
for i = 1:length(transitions)-1,
	target(transitions(i):transitions(i+1)) = repmat(i,transitions(i+1)-transitions(i)+1,1);
end

%  Reorder consolidated interval IDs
target(order) = target;
target = target';

consolidated = unique(intervals,'rows');

% Empty intervals belong to none
empty = diff(original,1,2) < 0;
target(empty) = NaN;
