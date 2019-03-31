%islmatrix - Test if parameter is a logical matrix (>= 2 columns).
%
%  USAGE
%
%    test = islmatrix(x)
%
%    x              parameter to test
%
%  SEE ALSO
%
%    See also islscalar, islvector, isdmatrix, isdvector, isdscalar, isivector, isiscalar,
%    isstring.
%

% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = islmatrix(x)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help islmatrix">islmatrix</a>'' for details).');
end

% Test: logical, two dimensions, two or more columns?
test = islogical(x) & length(size(x)) == 2 & size(x,2) >= 2;
