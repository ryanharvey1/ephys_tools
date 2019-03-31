function ts = findTime(D, ix)

% tsd/findTime  Returns timestamp at which specified index ix occurs
%
% ts = findTime(D, ix)
%
% INPUTS:
%       D - tsd object
%       ix - index of the timestamp you want to find
% OUTPUTS:
%       ts - timestamp at which index ix occurs
%
% ADR 1998, version L4.0, last modified '98 by ADR
% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

ts = D.t(ix);