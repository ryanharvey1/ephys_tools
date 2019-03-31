function ix = findAlignment(D, tstmp)

% ctsd/findAlignment  Returns index of specified timestamp in ctsd
%
% ix = findAlignment(D, tstmp)
%
% INPUTS:
% 	    D - ctsd object
%       tstmp - timestamp that you want to find the index of
% OUTPUTS:
%       ix - index of timestamp in D closest to tstmp
%
% ADR, version L4.0, last modified by ADR

% status: PROMOTED

ix = round((tstmp - D.t0)/D.dt) + 1;