function TS0 = Restrict(TS, t0, t1)

% ts/Restrict  Returns ts restricted to specified time interval(s)
%
% TS0 = Restrict(TS, t0, t1)
%
% INPUTS:
%       TS - ts object
%       t0, t1 - start and stop timestamps (possibly vectors) of desired intervals to include;
% OUTPUTS:
%       TSO - ts object which only includes times specified in inputs
%
% ADR, version L4.0, last modified by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

if size(t0) ~= size(t1)
   error('t0,t1 sizes do not match.');
end

f = [];
for iT = 1:length(t0)
   f = [f; find(TS.t >= t0(iT) & TS.t <= t1(iT))];
end
TS0 = ts(TS.t(f));
