function R = Range(tsa, tsflag)

% ctsd/Range  Returns time range covered by ctsd
%
% R = Range(tsa, tsflag)
%
% INPUTS:
%       tsa - ctsd object
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'sec0' returns time in sec counting from 0
%                if 'ms' returns time in ms
%                if 'all_ts', then range returns all timestamps in range
% OUTPUTS:
%       R - array of times covered by tsa
%
% ADR 1998, version L4.1, last modified 10/28/98 by ADR

% status: PROMOTED
% v4.1 28 oct 1998 flag no longer optional
% v4.2 23 Feb 2005 PL  flag optional again

if nargin == 1
    tsflag = 'ts';  % set default
end

R = StartTime(tsa):tsa.dt:EndTime(tsa);
R = R';
switch (tsflag)
    case 'sec'
        R = R/10000;
    case 'sec0'
        R = (R - StartTime(tsa))/10000;
    case 'ts'
    case 'ms'
        R = R/10;
    case 'all_ts'
        R = repmat(R, [1 tsa.dt]);
    otherwise
        error('Range called with invalid tsflag.');
end


