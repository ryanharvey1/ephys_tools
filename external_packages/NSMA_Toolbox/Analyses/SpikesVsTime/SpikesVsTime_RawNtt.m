function SpikesVsTime_RawNtt(tetrode_no, binsize, epochs)
% SpikesVsTime_RawNtt  Creates graph of spike rate vs. time using raw uncut
% Ntt file.
%
% SpikesVsTime_RawNtt(tetrode_no, binsize, epochs)
%
% tetrode_no = number of tetrode you wish to analyze (e.g., enter 10 for
%              Sc10.Ntt).  Currently code only recognizes files in the 
%              new neuralynx format Sc*.Ntt.
% binsize = size of bin (in msec) used for averaging data (optional)
% epochs = strucure containing names and time intervals.  See documention
%          for DrawEpochLines for details, or look at example below.  If
%          entered as an argument, intevals should be in microseconds.
%          For backwards compatibility you can also save your epoch
%          structure to a file called "epochs.mat" in the same directory as
%          your data.  In that case, intervals should be in old timestamps
%          (microseconds/100).
%
%
% Produces "Spikes Rate vs. Time" graph for a tetrode.  All spike events in
% the Sc*.Ntt file are included in this rate.
% Run from folder containing raw .Ntt Sc data files.
% Pass "epochs" variable as an argument or make sure 'epochs.mat' is saved 
% in same folder as your spike data.
% The function 'DrawEpochLines' should be in your path.
% Example:
%        epochs.intervals = {[3217254467 4626530863]/100, [5464343236 6543322119]/100};
%        epochs.names = {'maze1', 'maze2'};
%        save epochs.mat epochs 
%
%        SpikeVsTime_RawNtt(3);  % this will make a graph of Sc3.Ntt


% DE 3/30/07


if nargin<2
    binsize = 100;  % in ms
end;

% if user passes in epochs as a command line arguement, convert timestamps
% from microseconds into old-style 100 microsecond units needed by
% DrawEpochLines
if exist('epochs', 'var')
    for epi = 1:length(epochs.intervals)
        epochs.intervals{epi} = epochs.intervals{epi}/100;
    end;
end;

if nargin<3
    % get epochs from file
    if exist('epochs.mat')
        load epochs
    end %if exist('epochs.mat')
end;

% look in present working directory for files
spikefiles_temp = FindFiles(['Sc' num2str(tetrode_no) '.ntt']);
if isempty(spikefiles_temp)
    spikefiles_temp = FindFiles(['TT' num2str(tetrode_no, '%02g') '*.ntt']);
else   
    error(['Could not find spike file: ' 'Sc' num2str(tetrode_no) '.ntt']);
end;
        
% set up the fields used by Nlx2MatSpike
FieldSelection(1) = 1;
FieldSelection(2) = 0;
FieldSelection(3) = 0;
FieldSelection(4) = 0;
FieldSelection(5) = 0;
ExtractHeader = 1; % yes, do
ExtractionMode = 1;  % Extract All Spike Times

% note we subtract one from the look-up indices because Nlx2MatSpike
% uses C-style indices (starting from 0)
[ts, NLX_header] = ...
  Nlx2MatSpike(spikefiles_temp{1}, FieldSelection,ExtractHeader, ExtractionMode);
      
min_ts = min(ts);
max_ts = max(ts);

bins = min_ts:binsize*1e3:max_ts;
fr = hist(ts, bins);
fr_hz = fr/(binsize/1e3);

plot(bins/(60*1e6), fr_hz, 'k');
ax = axis;
hold on; DrawEpochLines(ax(4),0,epochs);

xlabel('time (min)');
ylabel('firing rate (Hz)');
title(['firing rate vs time, Tetrode ' num2str(tetrode_no)]); 
hold off;

return;
  