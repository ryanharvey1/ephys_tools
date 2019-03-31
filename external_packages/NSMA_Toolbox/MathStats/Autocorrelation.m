function Autocorrelation(TS)

% Autocorrelation(TS)
%
% INPUTS
%    TS = ts object
%
% Plots autocorrelation function

% ADR 1998
% version L4.0
% Status PROMOTED
ts = Data(TS);
[acorr,lags] = AutoCorr(ts, 10, 100);  % (timestamps, binwidth[msec], nbins)

bar(lags/1000, acorr);			% show acorr
H = findobj(gca, 'Type', 'patch');
set(H, 'facecolor', [0 0 0])
set(gca, 'YLim', [0 1.1*max(acorr)]);
title('autocorrelation')
