function bn = BinnedIntervals(start_and_end_times, dt)
%  Divide a list of intervals into uniform bins
%   bn = BinnedIntervals(start_and_end_times, dt)
%
% INPUT: [n,2] two-column array of start and end times of a list of intervals
%        and a bin size dt (in same units as start/end times).
%
% OUTPUT:
%        A two column array of start/end times of each bin where the intervals in
%        start_and_end_times are divided into bins of size dt. The first bin starts eaxctly at
%        the start time of each interval. The last bin ends WITHIN or at the end of each interval.
%        It never goes outside the endtime! Any bins that partially go beyond the end times are eliminated.
%        End-times have a slight amount (a matlab eps) subtracted so that they are numerically 
%        not equivalent to the start time.
% Original code by SC (called binned.m);
% slightly modified and renamed to BinnedIntervals.m by PL 2008
bn = [];
for ii = 1:size(start_and_end_times,1)
    s = start_and_end_times(ii,1):dt:start_and_end_times(ii,2);
    e = s + dt - eps;
    if e(end) > start_and_end_times(ii,2)
        e(end) = [];
        s(end) = [];
    end
    bn = [bn; s(:) e(:)];
end
