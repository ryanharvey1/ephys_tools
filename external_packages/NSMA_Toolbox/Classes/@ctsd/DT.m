function dt = DT(tsa)

% ctsd/DT  Returns DT (timestep between sequential entries) from ctsd
%
% dt = DT(tsa, tsflag)
%	
% INPUTS:
%       tsa - ctsd object
%       tsflag - if 'ts' returns time in timestamps (default),
%                if 'sec' returns time in sec
%                if 'ms' returns time in ms
% OUTPUTS:
%       dt - the tsa.dt (timestep between sequential data entries)
%
% ADR, version L4.0, last modified by ADR

% status: PROMOTED

dt = tsa.dt;

if nargin == 2
   switch tsflag
   case 'sec'
      dt = dt/10000;
   case 'ms'
      dt = dt/10;
   case 'ts'
   otherwise
      error('Unknown tsflag.');
   end
end
