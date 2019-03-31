function T1 = EndTime(tsa, tsflag)

% tsd/EndTime  Returns final timestamp in tsd 
%
% T1 = EndTime(tsa, tsflag)
%
% INPUTS:
%       tsa - tsd object 
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'ms' returns time in ms
% OUTPUTS:
%       T1 - final timestamp in tsa
%
% ADR 1998, version L4.0, last modified '98 by ADR
% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

T1 = max(tsa.t);
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