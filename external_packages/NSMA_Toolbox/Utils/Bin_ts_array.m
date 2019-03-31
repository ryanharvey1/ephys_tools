function [O,R] = Bin_ts_array(ts_array, start_end_times);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [O,R] = Bin_ts_array(ts_array, start_end_times);
%
% A function to bin a cell array of timestamps into a matrix. Unlike other 
% binning programs, the bin start and end times are completely up to the caller.
% 
% INPUT: a cell array of timestamps(as vectors) or ts objects
%        a n x 2 matrix of start and end time to bin the timestamps. Bins are
%         made from col1 up to but NOT including col 2. Intervals CAN overlap.
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
% cowen 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncells = length(ts_array);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to vectors if necessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(ts_array{1},'ts')
    for ii = 1:ncells
        ts_array{ii} = Data(ts_array{ii});
    end
end
if length(start_end_times) == 1
    % user passed a dt_msec: Create your own start_and_end_times.
    dt = start_end_times;
    mn = inf;
    mx = 0;
    
    for ii = 1:ncells
        if ~isempty(ts_array{ii})
            mn = min([ts_array{ii}(1) mn]);
            mx = max([ts_array{ii}(end) mx]);
        end
    end    
    if (~isinf(mn) & mx~=0)
        tt = mn:dt:mx;
        start_end_times = [tt(:) tt(:)+dt];
    else
        disp('Could not find any spikes')
        O = [];
        return
    end
    %disp('Created range.')
end

nbins = size(start_end_times,1);
O = zeros(nbins,ncells);

for ii = 1:ncells
    O(:,ii) = bin_times_by_intervals(ts_array{ii},start_end_times);
    %fprintf('.');
end

if nargout == 2
    R = start_end_times;
end