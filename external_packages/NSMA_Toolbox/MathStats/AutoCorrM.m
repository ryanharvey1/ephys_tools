function [C,B] = AutoCorrM(t1, binsize,nbins)
%function [C,B] = AutoCorrM(t1, binsize,nbins)
%
%  - called the same as Peter's AutoCorr.c except uses CrossCorr.c.
%  This function was created because AutoCorr was producing seg faults and I could not track
%  down the cause. (cowen). The faults were intermittent.
%   
%  INPUT: t1: a time series to auto correlate in 1/10000 sec
%                  (assumed to be sorted) 
%          binsize: the binsize for the auto corr histogram in msec
%          nbins: the number of bins
%  OUTPUT: C the auto correlation histogram
%           B (optional) a vector with the times corresponding to the bin centers
%  
%
% cowen
[C,B] = CrossCorr(t1,t1, binsize, nbins*2);
% Get rid of the large center bin.
mid = ceil(length(C)/2)+1;
C = C(mid:end);
% x is the center of each bin.
B = B(mid:end);% + binsize/2 ;
