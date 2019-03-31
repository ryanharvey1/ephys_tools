function [bi,cor,lag]=burstindex(spk,t_bin)
% Burst Index: (Coletta et al. 2018)
% Burst Index was defined as the sum of spikes with an interspike interval<
% 6ms, divided by the number of spikes

% Input: spk: timestamps in seconds
% Output: bi: Burst Index

% Ryan Harvey

max_lag = 0.5;
if ~exist('t_bin','var')
    t_bin=0.002; % 2ms
end
% Acor - taken from intrinsic frequency 2
if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
    max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
end
[cor, lag] = CrossCorr(spk, spk, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', 'count');

% find 0ms
zero=find(lag==0);
% find < 6ms
sixms=find(lag<0.006);
sixms=sixms(end);

bi=sum(cor(zero:sixms,1))/sum(cor(sixms+1:end));
end

