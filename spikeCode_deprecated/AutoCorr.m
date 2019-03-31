function [C, B] = AutoCorr(t1, binsize, nbins)
% [C, B] = AutoCorr(t1, t2, binsize, nbins)
%
% factorial Auto Correlation of a spike train time series.
%
% The 'factorial' version of Auto Correlation does not double-count
% the spike at lag 0 (i.e. it only counts pairs of DIFFERENT spikes, but
% ignores pairings of each spike with itself). This results in a smaller 
% peak in the 0-lag bin compared to calculating the autocorr via the
% crosscorr(t1,t1, binsize, nbins) function.
%
% INPUTS
% t1: array containing sorted time series
% binsize: size of the bin for the cross correlation histogram in msec
% nbins: number of bins in the histogram
% OUTPUTS
% C: the auto correlation histogram (only the right half - including the 0 lag - of the symetric
%                                     function is computed)
% B: a vector with the times corresponding to the bins (optional)

% batta 1999, lipa 2000
% MEX file
% status: beta/* AutoCorr
