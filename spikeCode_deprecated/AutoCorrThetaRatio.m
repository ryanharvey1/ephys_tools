% function AutoCorrThetaRatio(path, binsize, nbins)
%
% Input: 
%       - binesize: msec/bin (default is 4)
%       - nbins: number of bins (default is 250)
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
path = '/Users/bjclark/Documents/MATLAB/test_data/spike_data/2011-04-15_08-58-00/';

% find spike timestamp data
cell_list = FindFiles('TT*SS*.txt', 'StartingDirectory', path);

% loop through spike data and generate autocorr for each cell
for i = 1:length(cell_list);
    rawData = importdata(cell_list{i});
    [h, t] = AutoCorr(rawData, binsize, nbins);
    figure (i), plot(t, h, 'k'), axis tight;
    title('Autocorrelation');
    ylabel('rate')
    xlabel('msec (4 msec binsize)')
    xlim([0 1000]);
    box off
    [filepath, filename] = fileparts(cell_list{i});
    save([filepath filesep filename '.mat']);
    print(figure (i), '-djpeg', filename);
    close(figure (i))
    clear rawData h t
end



