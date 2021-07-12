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
ts = ts(1:length(vel),1); % make sure time is equal to velocity
tsidx=discretize(ts,ts(1):binsize:ts(end));
% mean rows together which correspond to each bin
n=min(tsidx):max(tsidx);
binnedvel=arrayfun(@(i) nanmean(vel(tsidx==i)),n)';
end