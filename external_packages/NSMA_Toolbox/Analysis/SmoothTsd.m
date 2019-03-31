function S = SmoothTsd(V, n)

% function S = SmoothTsd(V, n)
%
% Takes a tsd object and smooths the data according to an n-point
% symmetric hamming window.  Returns smoothed tsd.
%
% inputs:   V = tsd object
%           n = # of points (frames) that defines hamming window
%
% outputs:  S = tsd with smoothed data
%
% PL


n = floor(n/2) * 2;  % make # points even

hh = hamming(n);  % create n-point symmetric hamming window in a column vector
hh = hh/sum(hh);  % normalize to sum to 1 (make a true density function)

v = Data(V);
t = Range(V, 'ts');

v = conv(v, hh);  % take moving average at each data point using hamming window (smooths data)

v = v(n/2:end-n/2);     % chop off extra values added to beginning and end of data in
                        % conv process

S = tsd(t, v);  % outputs tsd with smoothed data
