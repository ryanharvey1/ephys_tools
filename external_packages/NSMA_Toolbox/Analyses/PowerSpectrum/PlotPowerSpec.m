function PlotPowerSpec(freqs, PowerAtFreqs, epochs, lo_freq, hi_freq, log_xaxis)

% PlotPowerSpec  Plots power spectrum by epoch - used internally within PowerSpec
%
% PlotPowerSpec(freqs, PowerAtFreqs, epochs, lo_freq, hi_freq, log_xaxis)
%
% INPUTS:
%   freqs - cell array with 3 cells corresponding to the 3 epochs (S1,M,S2); 
%            each cell contains a list of frequencies (Hz) over which the power spectrum is plotted
%   PowerAtFreqs - cell array w/ cells corresponding to epochs, as above.  Each cell contains list 
%                   of powers corresponding to list of frequencies for that epoch
%   epochs - a struct object with fields:
%       epochs.names = cell array of strings of espoch names 
%       epochs.intervals = cell array of 1x2 array with [start_ts  end_ts] (start
%           and end timestamps of each epoch) -- elements in cell array correspond
%           to epochs, listed in same sequence as in epochs.names
%   lo_freq, hi_freq - range of frequencies over which to plot power spectrum
%   log_xaxis - (optional input) enter the string 'log' if you want frequencies on
%       the x-axis plotted on a log_10 scale
% OUTPUTS:
%   (none)
%
% MRN 7/03, last modified 5/04 by MN


sc = [0 0 1 1];
fh = figure(1);
set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
orient tall;

min_y = min(min(PowerAtFreqs));
max_y = max(max(PowerAtFreqs));

if isempty(lo_freq) & isempty(hi_freq)
    lo_freq = min(min(freqs));
    hi_freq = max(max(freqs));
end %if isempty

localpower = PowerAtFreqs(find(freqs >= lo_freq & freqs <= hi_freq));
min_y = min(localpower);
max_y = max(localpower);

if strcmp(log_xaxis, 'log')
    xticks = [1:10 20 30 40 50 60 100 300];
    plot(log10(freqs),PowerAtFreqs);
    set(gca, 'XTick', log10(xticks), 'XTickLabel', xticks, 'YTick', [], 'YTickLabel', []); 
    min_x = log10(lo_freq);
    max_x = log10(hi_freq);
else
    plot(freqs,PowerAtFreqs);
    set(gca, 'YTick', [], 'YTickLabel', []);
    min_x = lo_freq;
    max_x = hi_freq;
end %if strcmp

legend(epochs.names);
axis([min_x max_x min_y max_y]);
xlabel('Frequency (Hz)');
ylabel('Power (cross-spectral density)');

   
% Save as Matlab fig file
saveas(fh,'PowerSpec.fig');
    
close all;
