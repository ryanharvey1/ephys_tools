function [M, SM, B, R] = PETH(X, Y, binsize_msec, nbins)

% PETH  Calculates Peri-Event-Time-Histogram
%
% [M, SM, B, R] = PETH(X, Y, binsize_msec, nbins)  
%
% INPUTS: 
%   X - array of timestamps (0.1 msec units) to be histogrammed
%   Y - array of trigger-event timestamps around which the histograms will be centered
%   binsize_msec - binsize in msec  
%   nbins - total number of bins centered around zero so that the histogram range is
%               [-(nbins*binsize_msec)/2, (nbins*binsize_msec)/2]
%               e.g., binsize_msec = 40;  nbins = 101,  <==>  (-2, +2) secs
%                     binsize_msec = 6;   nbins = 101,  <==>  (-300, +300) msec
%  OUTPUTS:
%   M - histogram values  (nbin x 1 vector)
%   SM - error of bin contents  
%   B - bins in msec
%   R - individual histogram contents around each event before averaging (nbin x length(Y) array) 
%
% PL 2001, last modified '02 by MN


switch nargin 
case 2
    %binsize = 40 * 10;   % {-2,+2) secs
    %binsize = 100 * 10;   % (-5, +5) sec
    binsize = 6 * 10;     % (-300, +300) msec
    nbins = 100 + 1;
case 3
    nbins = 100 + 1;
    binsize = binsize_msec*10;
case 4
    binsize = binsize_msec*10;
end

M = zeros(nbins, 1);
SM = zeros(nbins, 1);
B = -(nbins - 1) * binsize / 2: binsize:  (nbins - 1) * binsize / 2;
B = B / 10;
R = zeros(nbins, length(Y));

for i = 1: length(Y)
  edges = Y(i)- nbins * binsize / 2: binsize : Y(i) + nbins * binsize / 2;
  Q = histc(X, edges);
  Q = Q*10000/binsize;  % convert to spikes/sec
  M = M + Q(1:end-1);
  R(:, i) = Q(1:end-1);
  SM = SM + Q(1:end-1) .^2;
end

% standard error of mean spikes/sec
SM = sqrt( (SM-(M.^2)/length(Y)) / (length(Y)-1) / length(Y) );

% mean spikes/sec
M = M/length(Y);
