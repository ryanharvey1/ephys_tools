function T1 = StartTime(tsa, tsflag)

% ctsd/StartTime  Returns first timestamp in ctsd
%
% T1 = StartTime(tsa, tsflag)
%
% INPUTS:
%       tsa - ctsd object
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'ms' returns time in ms
% OUTPUTS:
%       T1 - first timestamp in tsa
%
% ADR, version L4.0, last modified by ADR

% status: PROMOTED

T1 = tsa.t0;

if nargin == 2
   switch tsflag
   case 'sec'
      T1 = T1/10000;
   case 'ms'
      T1 = T1/10;
   case 'ts'
   otherwise
      error('Unknown tsflag.');
   end
end
