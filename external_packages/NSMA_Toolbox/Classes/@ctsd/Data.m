function v = Data(tsa, ix)

% ctsd/Data  Retrieves data from ctsd
%
% v = Data(tsa)
% v = Data(tsa, ix)
%
% INPUTS:
%       tsa - ctsd object
%       ix - alignment list (timestamps)
% OUTPUTS:
%       v - the tsa.Data
%
%   if called with alignment list, returns those tsa.Data(ix)
%   if called without, returns complete tsa.Data
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status: PROMOTED

switch nargin
case 2
   f = findAlignment(tsa, ix);
   v = SelectAlongFirstDimension(tsa.data,f);
  case 1
    v = tsa.data;
  otherwise
    error('Unknown number of input arguments');
end
