function R = Restrict(D, t0, t1)

% ctsd/Restrict  Returns ctsd restricted to specified times
%
% R = Restrict(D, t0, t1)   [syntax 1]
% R = Restrict(D, t)        [syntax 2]
%
% INPUTS:
%       D - ctsd object
%       t0, t1 - [syntax 1] start and stop timestamps (possibly vectors) of desired intervals to include;
%       t - [syntax 2] list of timestamps to include;
% OUTPUTS:
%       R - ctsd or tsd object which only includes times specified in inputs
%
% ADR 1998, version L5.0, last modified 1/19/98 by ADR

% status: PROMOTED
% v4.1 29 oct 1998 now can handle nargin=2
% v5.0 30 oct 1998 time dimension is always 1st dimension
% v5.1 19 jan 1998 now can handle t0 and t1 as arrays

switch nargin
case 2                             % R = Restrict(tsd, t)
   ix = findAlignment(D, t0);
   tlist = Range(D, 'ts');
   R = tsd(tlist(ix), SelectAlongFirstDimension(D.data, ix));
   
case 3                             % R = Restrict(tsd, t0, t1)
   if length(t0) ~= length(t1)
      error('t0 and t1 must be the same length.');
   end
   if length(t0) == 1
      R0 = max(1, findAlignment(D, t0));
      R1 = min(length(D.data), findAlignment(D, t1));
      
      R = ctsd(t0, D.dt, SelectAlongFirstDimension(D.data, R0:R1));
   else
      DATA = [];
      TIME = [];
      for it = 1:length(t0)
         R0 = max(1, findAlignment(D, t0(it)));
         R1 = min(length(D.data), findAlignment(D, t1(it))); 
         TIME = cat(2, TIME, findTime(D, R0):D.dt:findTime(D, R1));
         DATA = cat(1, DATA, SelectAlongFirstDimension(D.data, R0:R1));
      end
      R = tsd(TIME, DATA);
   end
   
otherwise
   error('Unknown number of input arguments.');
   
end % switch
