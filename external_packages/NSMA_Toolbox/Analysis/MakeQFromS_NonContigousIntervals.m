function [Q,R] = MakeQFromS_NonContigousIntervals(S, start_end_times,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Q,R] = MakeQFromS_NonContigousIntervals(S, start_end_times,dt);
%
% A function to bin a cell array of timestamps into a matrix. Unlike other 
% binning functions (e.g. MakeQFromS), the user can give a list of non-contiguous intervals (which may even overlap - but is not recomended)
% and the intervals are binned such that the first bin startsts at the beginning of each interval and the last bin always lies FULLY within each interval.
% 
% INPUT: S = a cell array of timestamps(as vectors) or ts objects
%        and start_end_times =  n x 2 matrix of start and end time to bin the timestamps. The 3rd argument dt = binsize in same units as start_end_times).
%        Bins are made from col1 up to but NOT including col 2. Intervals CAN overlap.
%         OR, you can just pass a single value. This is assumed to be the space between
%         timestamps in whatever units the timestamps are in. 
%         Results will be binned from the first to the last timstamp with this bin size.
%
% OUTPUT: a matrix with rows = to the number of start and end times and cols
%         equal to the number of elements in the ts_array.
%         optionally the start and end times of each bin in R.
% 
%  NOTE: MakeQFromS bins on the centers of the values returned from the Range command whereas
%        Bin_ts_array bins between the timestamps returned in R. (>=col1 and < col2). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original code by cowen 
% Modified by PL June 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncells = length(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to vectors if necessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(S{1},'ts')
    for ii = 1:ncells
        ts_array{ii} = Data(S{ii});
    end
end
bn = BinnedIntervals(start_end_times, dt);
nbins = size(bn,1);
O = zeros(nbins,ncells);

for ii = 1:ncells
    O(:,ii) = bin_times_by_intervals(ts_array{ii},bn);
end

if nargout == 2
    R = bn;
end

Q = tsd(bn(:,1),sparse(O));