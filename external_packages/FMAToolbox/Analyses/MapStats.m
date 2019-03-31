%MapStats - Compute statistics for a map Z = f(X,Y), or a curve Z = f(X).
%
%  Compute statistics on a continuous map, where one time-varying variable Z
%  is represented as a function of one or two time-varying variables X and Y.
%  The variable Z can either be a point process (typically, a list of spike
%  timestamps) or a continuous measure (e.g. the instantaneous velocity of
%  the animal, the spectral power of an LFP channel in a given frequency band,
%  the coherence between two oscillating LFP channels, etc.) Typical examples
%  of X and Y include spatial coordinates and angular directions.
%
%  USAGE
%
%    stats = MapStats(map,<options>)
%
%    map            map obtained using <a href="matlab:help Map">Map</a>
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'center'      search for the center of the field around given
%                   coordinates (to manually select one of multiple fields)
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.2)
%     'minSize'     fields smaller than this size are considered spurious
%                   and ignored (default = 100 for 2D, 10 for 1D)
%     'type'        'll' if X and Y are linear, 'cl' if X is circular and Y
%                   linear, 'lc' if X is linear and Y circular, or 'cc' if X
%                   and Y are circular - for 1D data, a single letter is used
%                   (default = 'll')
%    =========================================================================
%
%  OUTPUT
%
%    stats.x             abscissa of the maximum value (in bins)
%    stats.y             ordinate of the maximum value (in bins)
%    stats.peak          in-field maximum value
%    stats.mean          in-field mean value
%    stats.size          field size (in bins)
%    stats.field         field (1 = bin in field, 0 = bin not in field)
%    stats.fieldX        field x boundaries (in bins)
%    stats.fieldY        field y boundaries (in bins)
%    stats.specificity   spatial specificity (Skaggs et al., 1993)
%
%  SEE
%
%    See also Map, FiringMap, FiringCurve.

% Copyright (C) 2002-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function stats = MapStats(map,varargin)

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
end

if ~isfield(map,'z'),
	map.z = map.rate;
	map = rmfield(map,'rate');
end

% Default values
stats.peak = 0;
stats.mean = 0;
stats.size = 0;
stats.field = [];
center = [];
threshold = 0.2;
type = 'll';

nDims = sum(size(map.z)>=2);
if nDims == 2,
	minSize = 100;
else
	minSize = 10;
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).']);
	end
	switch(lower(varargin{i})),

		case 'center',
			center = floor((varargin{i+1}+1)/binSize);
			if ~isivector(center,'#2','>0'),
				error('Incorrect value for property ''center'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end

		case 'threshold',
			threshold = varargin{i+1};
			if ~isdscalar(threshold,'>=0','<=1'),
				error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end

		case 'minsize',
			minSize = varargin{i+1};
			if ~isdscalar(minSize,'>=0'),
				error('Incorrect value for property ''minSize'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end

		case 'type',
			type = lower(varargin{i+1});
			if (nDims == 1 && ~isstring(type,'c','l')) || (nDims == 2 && ~isstring(type,'cc','cl','lc','ll')),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).');
			end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MapStats">MapStats</a>'' for details).']);

  end
end

% Are X and/or Y circular?
circX = size(map.z,2) > 1 && strcmp(type(1),'c');
circY = size(map.z,1) > 1 && ((size(map.z,2) > 1 && strcmp(type(2),'c')) || strcmp(type(1),'c'));

% Default values
stats.x = NaN;
stats.y = NaN;
stats.field = [];
stats.size = 0;
stats.peak = 0;
stats.mean = 0;
stats.fieldX = NaN;
stats.fieldY = NaN;

if isempty(map.z), return; end

% Compute the spatial specificity of the map, based on the formula proposed by Skaggs et al. (1993).
T = sum(sum(map.time));
if T == 0,
  stats.specificity = 0;
else
  occupancy = map.time/(T+eps);
  m = sum(sum(map.count))/(sum(sum(map.time)+eps));
  if m == 0,
    stats.specificity = 0;
  else
    logArg = map.count/m;
    logArg(logArg <= 1) = 1;
    stats.specificity = sum(sum(map.count.*log2(logArg).*occupancy))/m;
  end
end

% Determine the field (i.e. the connex area where the value or rate is > threshold*peak).
% There are two ways to do this:
% 1) If the location of the center was provided (see 'center' property), the field is
%    found around the center.
% 2) Otherwise, the field is assumed to be around the bin with maximal value or rate
%    (when there are two fields, the one with the highest value or rate is selected unless
%    the other one is >20% bigger).

% Circular X and/or Y
if circX,
	map.z = [map.z map.z];
	map.time = [map.time map.time];
	map.x = [map.x map.x+map.x(end)];
end
if circY,
	map.z = [map.z ; map.z];
	map.time = [map.time ; map.time];
	map.y = [map.y map.y+map.y(end)];
end

if max(max(map.z)) == 0,
  stats.field = logical(zeros(size(map.z)));
  return;
end

nBinsX = max([1 length(map.x)]);	% minimum number of bins is 1
nBinsY = max([1 length(map.y)]);

% First clean the map from very small (spurious) regions of elevated values or rates
map.z = CleanMap(nBinsX,nBinsY,threshold,minSize,map.z);

if isempty(center),
	% In the absence of any information regarding the location of the center,
	% find two candidate fields
	z = map.z;
	for i = 1:2,
		% Find field and compute all parameters
		peak(i) = max(max(z));
		peakLocation = FindPeakLocation(z);
		x(i) = peakLocation(1);
		y(i) = peakLocation(2);
		field{i} = FindFiringField(x(i),y(i),nBinsX,nBinsY,threshold*peak,z);
		fieldSize(i) = sum(sum(field{i}));
		% Remove this field from the map for next iteration
		z(logical(field{i})) = 0;
	end
	% Choose between the two candidate fields
	if fieldSize(2) == 0,
		winner = 1;
	else
		%% peakRelativeDifference = abs(peak(2)-peak(1))/((peak(1)+peak(2))/2);
		sizeRelativeDifference = (fieldSize(2)-fieldSize(1))/((fieldSize(1)+fieldSize(2))/2);
		if sizeRelativeDifference > 0.2,
			winner = 2; % Choose field #2 (see below)
		else
			winner = 1; % Choose field #1 (see below)
		end
	end
else
	% The approximate location of the center was provided in order to select one of multiple fields
	% First pass: find approximate field (threshold is not accurate since in-field peak is unknown)
	z = map.z;
	approximatePeak = z(center(2),center(1)); % (x,y) means column x and row y, thus rate(y,x)
	field{1} = FindFiringField(center(1),center(2),nBinsX,nBinsY,approximatePeak,z);
	% Remove all other potential fields from the rate map for next iteration
	z(logical(~field{1})) = 0;
	% Second pass: accurate field
	peak(1) = max(max(z));
	peakLocation = FindPeakLocation(z);
	x(1) = peakLocation(1);
	y(1) = peakLocation(2);
	field{1} = FindFiringField(x(1),y(1),nBinsX,nBinsY,threshold*peak,z);
	fieldSize(1) = sum(sum(field{1}));
	winner = 1;
end

% Set stats
stats.x = x(winner);
stats.y = y(winner);
stats.field = logical(field{winner});
stats.size = fieldSize(winner);
stats.peak = peak(winner);
stats.mean = mean(mean(map.z(stats.field)));
i = find(max(stats.field,[],1));
if ~isempty(i),
	stats.fieldX = [i(1) i(end)];
end
i = find(max(stats.field,[],2));
if ~isempty(i),
	stats.fieldY = [i(1) i(end)];
end

% Circular X and/or Y
if circX,
	stats.field = stats.field(:,1:end/2)|stats.field(:,end/2+1:end);
	if stats.fieldX(2) > nBinsX/2,
		stats.fieldX(2) = stats.fieldX(2) - nBinsX/2;
	end
end
if circY,
	stats.field = stats.field(1:end/2,:)|stats.field(end/2+1:end,:);
	if stats.fieldY(2) > nBinsY/2,
		stats.fieldY(2) = stats.fieldY(2) - nBinsY/2;
	end
end


% ------------------------------- Helper functions -------------------------------


% CleanMap - repeatedly find and discard small regions of elevated values or rates

function map = CleanMap(M,N,threshold,minSize,map,n)

if nargin < 6,
	n = 1;
elseif n >= 50,
	% Maximum 50 recursive calls
	return
else
	n = n + 1;
end

peak = max(max(map));
peakLocation = FindPeakLocation(map);
field = FindFiringField(peakLocation(1),peakLocation(2),M,N,threshold*peak,map);
fieldSize = sum(sum(field));
disp(['CleanMap: field size = ' int2str(fieldSize)]);
if fieldSize > 0 & fieldSize < minSize,
	% Remove this region from rate map for next iteration
	map(logical(field)) = 0;
	map = CleanMap(M,N,threshold,minSize,map,n);
end

% FindPeakLocation - find the coordinates of the peak value or rate

function peakLocation = FindPeakLocation(map)

peak = max(max(map));
xy = (map == peak);
[y,x] = find(xy);
peakLocation = [x(1) y(1)];
