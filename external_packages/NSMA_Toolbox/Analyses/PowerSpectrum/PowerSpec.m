function PowerSpec(lo_freq, hi_freq, log_xaxis)

% PowerSpec  Plots power spectrum of eeg data by epoch
%
% PowerSpec(lo_freq, hi_freq, log_xaxis)
%
% INPUTS:
%   lo_freq, hi_freq - range of frequencies over which to plot power spectrum
%   log_xaxis - (optional input) enter the string 'log' if you want frequencies on
%       the x-axis plotted on a log_10 scale
% OUTPUTS:
%   (none)
%
% Note: this function, particularly use of mtcsd (multitaper cross-spectral density)
% is based on my possibly incomplete understanding of the code that Francesco 
% Battaglia gave me just before he left.  Others more familiar with power spectrum 
% calculations may want to look over the details more carefully.
%
% MN 5/04, from FB, last modified 5/04 by MN


% Get epochs
if exist('epochs.mat')
    load epochs;
end %if exist('epochs.mat')


% Get eeg data
eegfiles=findfiles('CSC*');
[d, fname, ext] = fileparts(eegfiles{1});
eeg_tsd = ReadCR_tsd([fname ext]);
[st, et, sFreq, numrec] = ReadCR_get_info([fname ext]);


% Calculate power spectra by epoch
freqs = [];
PowerAtFreqs = [];

% find nearest power of 2 just above sampling freq.
nFFT = 2; i=1;
while nFFT < sFreq
    nFFT = 2^i;
    i=i+1;
end %while

for i=1:length(epochs.intervals)
    x=data(restrict(eeg_tsd,epochs.intervals{i}(1),epochs.intervals{i}(2)));
    [yo,fo] = mtcsd(x,nFFT,sFreq,nFFT,nFFT/2,4,'linear',7);
    freqs(:,i) = fo;
    PowerAtFreqs(:,i) = log10(abs(yo));
end %for


% Plot
if nargin < 3
    log_xaxis = 0;
end %if
if nargin <1
    lo_freq = [];
    hi_freq = [];
end %if
    
PlotPowerSpec(freqs,PowerAtFreqs,epochs,lo_freq,hi_freq,log_xaxis);

