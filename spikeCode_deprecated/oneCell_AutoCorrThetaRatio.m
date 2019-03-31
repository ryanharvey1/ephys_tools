% function AutoCorrThetaRatio(path, binsize, nbins)
%
% Input: 
%       - binesize: msec/bin (default is 25)
%       - nbins: number of bins (default is 40)
%       - path: directory to folder with spike files
%
% Output:
%       - figure showing autocorrelation
%       - h: y-axis data
%       - t: x-axis
%
% Created by Ben C 2013; patched together from AutoCorr.m from MClust
%

% define binsize and number of bins
binsize = 4;
nbins = 250;

% define the path to spike sorted files
path = 'C:\Matlab2\EC_spikes';

% find spike timestamp data
cell_list = FindFiles('TT*SS*.txt', 'StartingDirectory', path);

% generate autocorr for cell
rawData = importdata(cell_list{1});
[h, t] = AutoCorr((rawData/100), binsize, nbins);
hi = h./max(h);
figure (1), plot(t, hi, 'k'), 
ylabel('Normalized Auto-Correlation')
xlabel('Lag (msec)')
xlim([0 500]);
box off


