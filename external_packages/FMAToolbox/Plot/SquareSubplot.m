%SquareSubplot - Layout subplots in a square arrangement.
%
%
%  USAGE
%
%    p = SquareSubplot(n,i)
%
%    n              total number of subplots
%    i              current subplot

% Copyright (C) 2008 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function p = SquareSubplot(n,i)

if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SquareSubplot">SquareSubplot</a>'' for details).');
end

M = round(sqrt(n));
N = ceil(n/M);
subplot(M,N,i);
