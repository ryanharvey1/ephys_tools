% * binning by interval
% * MEX file
% * 
% * cowen 2002
% * 
% * input: 1: Spike or other timestamped data -- n_spikes x 1, assumed to be sorted
% * input: 2: n_intervals x 2 matrix of intervals. The spike times will be binned into 
% *			 these intervals. Also assumed to be sorted. The spikes are counted from 
% the values in col 1 up to but NOT including the values in col 2.
% * output: A vector of length n_intervals that contains the spike counts to each element in the bin.
% 
%   * NOTE: ASSUMES TIMESTAMPS ARE SORTED AND THE INTERVALS ARE SORTED BY THE FIRST COLUMN.
