function R = Restrict(tsa, t0, t1)
% tsd/Restrict  Returns tsd restricted to specified times
%
% R = Restrict(tsa, t0, t1)   [syntax 1]
% R = Restrict(tsa, t)        [syntax 2]
%
% INPUTS:
%       tsa - tsd object
%       t0, t1 - [syntax 1] start and stop timestamps (possibly vectors) of desired intervals to include;
%       t - [syntax 2] list of timestamps to include;
% OUTPUTS:
%       R - tsd object which only includes times specified in inputs
%
% ADR, version L4.1, last modified 1/23/99 by ADR
% status: PROMOTED
% v4.1 29 oct 1998 now can handle nargin=2
% v4.2 23 jan 1999 now can handle t0 and t1 as arrays

switch nargin
case 2                             % R = Restrict(tsd, t)
   ix = findAlignment(tsa, t0);
case 3                             % R = Restrict(tsd, t0, t1)
   if length(t0) ~= length(t1)
      error('t0 and t1 must be same length')
   end
   ix = [];
   for it = 1:length(t0)
      f = find(tsa.t >= t0(it) & tsa.t <= t1(it));
      ix = cat(1, ix, findAlignment(tsa, tsa.t(f)));
   end
otherwise
   error('Unknown number of input arguments.');
end % switch
R = tsd(tsa.t(ix), SelectAlongFirstDimension(tsa.data, ix));
