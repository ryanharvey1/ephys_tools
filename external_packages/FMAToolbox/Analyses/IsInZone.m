%IsInZone - Test when the animal is in a given zone.
%
% Test when the animal is in a given zone of the experimental setup (typically the
% firing field of a place cell).
%
%  USAGE
%
%    status = IsInZone(positions,zone)
%
%    positions      position samples (normalized between 0 and 1)
%    zone           test area, obtained using e.g. <a href="matlab:help FiringMap">FiringMap</a> or <a href="matlab:help DefineZone">DefineZone</a>
%
%  OUTPUT
%
%    status         a logical vector indicating whether the animal was within
%                   (value 1) or outside (value 0) the zone at each timestamp
%
%  NOTE
%
%    In matrix 'zone', columns (resp. rows) correspond to x values (resp. y values)
%    for position samples.
%
%  SEE
%
%    See also FiringMap, MapStats, DefineZone.


% Copyright (C) 2004-2009 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function status = IsInZone(positions,zone)

% Check number of parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help IsInZone">IsInZone</a>'' for details).');
end

% Check parameter sizes
if size(positions,2) < 3,
	error('Parameter ''positions'' should have at least 3 columns (type ''help <a href="matlab:help IsInZone">IsInZone</a>'' for details).');
end
x = positions(:,2);
y = positions(:,3);
if max(x) > 1 || min(x) < 0,
	error('X coordinates should contain values in [0 1].');
end
if max(y) > 1 || min(y) < 0,
	error('Y coordinates should contain values in [0 1].');
end
% Discretize
X = size(zone,2);
Y = size(zone,1);
x = Bin(x,[0 1],X);
y = Bin(y,[0 1],Y);
% Now (x,y) are matrix indices in 'zone'. Convert them to linear indices
i = y+(x-1)*Y;
% Return timestamps for those indices that fall in the zone
status = zone(i);
