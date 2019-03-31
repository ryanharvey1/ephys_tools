%  [C,B,CErr] = CrossCorrWithStd(t1,t2,dx,nbins)
%
%  t1,t2 ... input timestamp vectors of cell 1 and 2 (in 0.1 msec units!!!) 
%            (must be sorted in ascending order!)
%  dx    ... lag binsize in msec!!!
%  nbins ... number of lag bins (should be odd - if nbins is even, nbins+1 C-values are calculated and returned) 
%            
%
%  C ...    cross-correlation or PETH of t2 around events given at t1 (t1 = events, t2 = snapshots) 
%  B ...    (optional) array of timestamps of lower-edges of the bins 
%  CErr ... (optional) a parallel vector to C holding the errors (standard deviations) of C
%              Errors are calculated as sqrt(n) where n = # of coincidences in each bin and divided by same
%              normalizing factor (nt1 * binsize / 10000) as the bins of C.
%  
%%  NOTES: 
%   1) Roles of t1 and t2: 
%       In Francesco's original version the CrossCorr function normalized such that it is
%       equivalent to a PETH (Peri-Event-Time-Histogram) of spike-train t2
%       around trigger-events given by spike-train t1. In this convention,
%       if spikes of t1 consistenly precede spikes of t2 by an interval -DT, C shows a peak to
%       the RIGHT of histogram center (with peak center at +DT). WARNING: This is
%       the opposite of the result obtained with matlabs bultin xcorr-function on binned
%       spiketrains!
%
%   2) Returned Bin Edges:
%       The bins are always calculated such that there is the same number of bins to
%       the left and right of 0. There is NO CENTER BIN spanning the 0-lag point. The
%       0-point ALWAYS lies exatly on a bin left-edge (e.g. if you enter nbins=200,
%       you get 201 bins returned and the bin 101 carries the counts of the interval [0,dx] to the right of the origin 
%       NOTE: to plot the resulting crosscorr histogram with matlabs bar function use:
%           >> bar(B,C,'histc'); 
%       to plot the histogram accurately. If you use the plot function use:
%           >> plot(B+dx/2,C);    
%       or
%           >> errorbar(B+dx/2,C,CErr);
%
%   3) VERY IMPORTANT:
%      t1 and t2 MUST be sorted in strictly ascending order! The alorithm produces
%      WILDLY WRONG results if even just a few spikes are out of sequence!!!
%
% * version 1.10
% *  added simple error calculation (PL & DE, June 2008)
% * version 1.01
% *  double added to prevent integer roundouff problem (cowen/battaglia)
% -----------------------------------*/

%% Original Mex File Docu:
% * cross correlations with error bars (+- 1 standard dev) in each bin
% * MEX file
% * 
% * batta 1999
% * 
% * input: t1, t2: two time series to cross correlate in 1/10000 sec
% *                (assumed to be sorted) 
% *        binsize: the binsize for the cross corr histogram in msec
% *        nbins: the number of bins
% * output: C the cross correlation histogram   
% *         B (optional) a vector with the times corresponding to the bins
% *         Cerr (optional) a parallel vector to C holding the errors (standard deviations) of C
% *              Errors are calculated as sqrt(n) where n = # of coincidences in each bin and divided by same
% *              normalizing factor (nt1 * binsize / 10000) as the bins of C.
% *
% * version 1.10
% *  added simple error calculation (PL & DE, June 2008)
% * version 1.01
% *  double added to prevent integer roundouff problem (cowen/battaglia)
% -----------------------------------*/
% % 