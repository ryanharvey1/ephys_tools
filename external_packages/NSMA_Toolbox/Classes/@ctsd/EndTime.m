function T1 = EndTime(tsa, tsflag)

% ctsd/EndTime  Returns final timestamp in ctsd 
%
% T1 = EndTime(tsa, tsflag)
%
% INPUTS:
%       tsa - ctsd object 
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'ms' returns time in ms
% OUTPUTS:
%       T1 - final timestamp in tsa
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status: PROMOTED

T1 = tsa.t0 + tsa.dt * (length(tsa.data)-1);

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
