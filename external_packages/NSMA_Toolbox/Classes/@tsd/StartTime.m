function T1 = StartTime(tsa, tsflag)

% tsd/StartTime  Returns first timestamp in tsd
%
% T1 = StartTime(tsa, tsflag)
%
% INPUTS:
%       tsa - tsd object
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'ms' returns time in ms
% OUTPUTS:
%       T1 - first timestamp in tsa
%
% ADR, version L4.0, last modified by ADR
% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

T1 = min(tsa.t);

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