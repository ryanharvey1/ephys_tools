% binsearch  (MEX FILE) Does a binary search to find index of specified value in a sorted list of data
%
% ix = binsearch_floor(data, key)
%
% INPUTS:
%       data - list of data.  Function assumes data is sorted in ascending order.
%       key - data value you want to find index of
% OUTPUTS:
%       ix - index of greatest value in data that is <= key  (i.e. ...<= data(ix-1) <= data(ix) <= key)
%               (if you need ix to be the value of the data value CLOSEST to key 
%               use the mex-file binsearch instead!)
%
% NOTE: if a key is outside the range of the data [data(1), data(end)] then ix = 1 or
%       ix = length(data), depending on which side of the data-range the key is
%       located. Test for (?key inside data range?) if you need different edge behaviour!    
%
% PL 2002 based on ADR 1998, version L4.0, last modified '98 by ADR
% PL 02/2009   corrected edge behaviour

