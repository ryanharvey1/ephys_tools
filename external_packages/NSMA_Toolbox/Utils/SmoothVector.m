function Z = SmoothVector(X,sigma)
% Smooth a vector X of numbers (row or column vector of lenght N) with a gaussian
% of width sigma (scalar) given in units of bins (e.g. a gaussian of half-width
% sigma = 3 would smooth over 3 bins to the left and right of the center
% bin)
% Z is a smoothed vector of same shape and lenght as X (row or column vector).
% Smoothed values at the left and right edge of Z imply that X is extended
% with zeros to the left and right ('zero padding').
%
% PL June 15,2007

% make gaussian filter of approriate length
N = length(X);
L = ceil(3*sigma);  % gaussian filter half length; cut the gaussian at 3 sigma
y = (-L:L)';  % filter length is 2L+1; column vector
gaussfilt = exp(-(y.^2/sigma^2));
gaussfilt = gaussfilt/sum(gaussfilt);  % gaussian filter normalized to sum = 1; column vector

% convolve signal X with gaussian filter
Z = conv(X,gaussfilt);   % Z now has length L + N + L
Z = Z(L+1:end-L);        % restrict smoothed vector Z to central N samples
 
