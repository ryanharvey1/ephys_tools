function t = StartTime(TS, tsflag)

% ts/StartTime  Returns first timestamp in ts
%
% t = StartTime(TS, tsflag)
%
% INPUTS:
%       TS - ts object
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'ms' returns time in ms
% OUTPUTS:
%       t - first timestamp in TS
%
% ADR 1998, version L4.1, last modified '98 by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m


if isempty(TS.t) 
   t = NaN;
else
   t = min(TS.t);

   if nargin == 2
       
      switch tsflag
      case 'sec'
         t = t/10000;
      case 'ms'
         t = t/10;
      case 'ts'
      otherwise
         error('Unknown tsflag.');
      end
      
   end
   
end