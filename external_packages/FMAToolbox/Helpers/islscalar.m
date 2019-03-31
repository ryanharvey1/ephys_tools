%islscalar - Test if parameter is a logical scalar.
%
%  USAGE
%
%    test = islscalar(x)
%
%    x              parameter to test
%
%  SEE ALSO
%
%    See also islvector, islmatrix, isdmatrix, isdvector, isdscalar, isimatrix, isivector,
%    isstring.
%

% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = islscalar(x)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help islscalar">islscalar</a>'' for details).');
end

% Test: double, scalar
test = islogical(x) & isscalar(x);
