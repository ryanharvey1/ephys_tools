function binnedvel=binned_vel(ts,vel,binsize)
% binned_vel: bins instant velocity 
%   Input:
%           ts: time stamps in seconds
%           vel: instant velocity
%           binsize: time bins in seconds 
%   Output:
%           binned_vel: binned velocity   
%
% Ryan Harvey 2019

% find index of rows for each bin
tsidx=discretize(ts,ts(1):binsize:ts(end));

% mean rows together which correspond to each bin
n=min(tsidx):max(tsidx);
binnedvel=arrayfun(@(i) nanmean(vel(tsidx==i)),n)';
end