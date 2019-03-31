% binsearch  (MEX FILE) Does a binary search to find index of specified value in a sorted list of data
%
% ix = binsearch(data, key)
%
% INPUTS:
%       data - list of data.  Function assumes data is sorted in ascending order.
%       key - data value you want to find index of
% OUTPUTS:
%       ix - index of value in data that is closest to key
%               (if you need ix to be the value of the immediate predecessor 
%               data(ix) <= key, use the mex-file binsearch_floor instead!)
%
% ADR 1998, version L4.0, last modified '98 by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in Contents.m