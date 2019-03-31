function bn = binned(start_and_end_times, dt);
%
%  bn = binned(start_and_end_times, dt)
%
% INPUT: n,2 matrix of start and end times
%        a bin size.
% OUTPUT:
%        start and end times of each bin where the intervals in
%        start_and_end_times are divided into bins of size dt. Any
%        bins that go beyond the end times are eliminated.
%        end times have a slight amount subtracted so that they are 
%        not equivalent to the start time.
%

bn = [];
for ii = 1:size(start_and_end_times,1)
    s = start_and_end_times(ii,1):dt:start_and_end_times(ii,2);
    e = s + dt - .0000001;
    if isempty(e)
        continue
    end
    if e(end) > start_and_end_times(ii,2)
        e(end) = [];
        s(end) = [];
    end
    bn = [bn; s(:) e(:)];
end
