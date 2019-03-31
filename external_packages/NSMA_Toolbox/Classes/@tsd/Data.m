function v = Data(tsa, ix)

% tsd/Data  Retrieves data from tsd
%
%   d = Data(tsa)
%   d = Data(tsa, ix)
%
% INPUTS:
%       tsa - tsd object
%       ix - alignment list (timestamps)
% OUTPUTS:
%       v - the tsa.Data
%
%   if called with alignment list, returns those tsa.Data(ix)
%   if called without, returns complete tsa.Data
%
% ADR 1998, version L4.1, last modified '98 by ADR
% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

switch nargin
case 2
   f = findAlignment(tsa, ix);
   v = SelectAlongFirstDimension(tsa.data,f);
case 1
   v = tsa.data;
otherwise
   error('Unknown number of input arguments');
end