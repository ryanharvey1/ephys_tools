function ts = findTime(D, ix)

% ctsd/findTime  Returns timestamp at which specified index ix occurs
%
% ts = findTime(D, ix)
%
% INPUTS:
%       D - ctsd object
%       ix - index of the timestamp you want to find
% OUTPUTS:
%       ts - timestamp at which index ix occurs
%
% ADR, version L4.0, last modified by ADR

% status: PROMOTED

ts = D.t0 + ix * D.dt;
