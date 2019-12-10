function X = wrap(Y, bound)
% X = WRAP(Y)
%   Maps Y values onto circular space (-PI <= X < PI). 
%   X = WRAP(Y, BOUND) specifies alternative bounds (default is PI).
%
%   --> www.paulbays.com

if nargin<2, bound = pi; end

X = mod(Y + bound, bound*2) - bound;
