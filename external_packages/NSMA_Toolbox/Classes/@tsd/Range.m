function R = Range(tsa, tsflag)

% tsd/Range  Returns time range covered by tsd
%
% R = Range(tsa, tsflag)
%
% INPUTS:
%       tsa - tsd object
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'sec0' returns time in sec counting from 0
%                if 'ms' returns time in ms
%                if 'all_ts', then range returns all timestamps in range
% OUTPUTS:
%       R - array of times covered by tsa
%
% ADR, version L4.1, last modified by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

if nargin == 1
    tsflag = 'ts';      % set default
end


switch (tsflag)
    case 'sec'
        R = tsa.t/10000;
    case 'sec0'
        R = (tsa.t - StartTime(tsa))/10000;
    case 'ts'
        R = tsa.t;
    case 'ms'
        R = tsa.t/10;
    case 'all_ts'
        R = StartTime(tsa.t):EndTime(tsa.t);
    otherwise
        error('Range called with invalid tsflag.');
end
