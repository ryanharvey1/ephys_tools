% oneCell_WaveFormCharacteristics
%
% Input: 
%       - AvgWaveform output from Nlx
%
% Output:
%       - Dur = the duration of waveforms on each wire
%       - meanDur = mean spike duration across all wires
%       - maxPeak
%
% Created by Ben C Jan 2014
%


% identify path to data file
path = '/Users/bjclark/Documents/MATLAB/test_data/spike_data/waveforms';
ReadData = FindFiles('*.csv', 'StartingDirectory', path);

% load data (first row is the timestamp; the other four rows are avg
% waveforms for each wire in order (1 to 4)
wvf = load(ReadData{1});

% convert to usec and uV
wvf = wvf*1000;

% wire 1 duration
w1_max = max(wvf(2,:));
w1_maxI = find(wvf(2,:) == w1_max);
w1max = wvf(1,w1_maxI);
w1max = w1max(1,1);

w1_min = min(wvf(2,:));
w1_minI = find(wvf(2,:) == w1_min);
w1min = wvf(1,w1_minI);
w1min = w1min(1,1);

w1_MinMax = [w1min w1max];
w1_Dur = max(w1_MinMax) - min (w1_MinMax);

% wire 2 duration
w2_max = max(wvf(3,:));
w2_maxI = find(wvf(3,:) == w2_max);
w2max = wvf(1,w2_maxI);
w2max = w2max(1,1);

w2_min = min(wvf(3,:));
w2_minI = find(wvf(3,:) == w2_min);
w2min = wvf(1,w2_minI);
w2min = w2min(1,1);

w2_MinMax = [w2min w2max];
w2_Dur = max(w2_MinMax) - min (w2_MinMax);

% wire 3 duration
w3_max = max(wvf(4,:));
w3_maxI = find(wvf(4,:) == w3_max);
w3max = wvf(1,w3_maxI);
w3max = w3max(1,1);

w3_min = min(wvf(4,:));
w3_minI = find(wvf(4,:) == w3_min);
w3min = wvf(1,w3_minI);
w3min = w3min(1,1);

w3_MinMax = [w3min w3max];
w3_Dur = max(w3_MinMax) - min (w3_MinMax);

% wire 4 duration
w4_max = max(wvf(5,:));
w4_maxI = find(wvf(5,:) == w4_max);
w4max = wvf(1,w4_maxI);
w4max = w4max(1,1);

w4_min = min(wvf(5,:));
w4_minI = find(wvf(5,:) == w4_min);
w4min = wvf(1,w4_minI);
w4min = w4min(1,1);

w4_MinMax = [w4min w4max];
w4_Dur = max(w4_MinMax) - min (w4_MinMax);

Dur = [w1_Dur w2_Dur w3_Dur w4_Dur];
meanDur = mean(Dur);

% find waveform peak amplitude
peak = w1_max;
peak(:,2) = w2_max;
peak(:,3) = w3_max;
peak(:,4) = w4_max;
maxPeak = max(peak);
PlotmaxPeak = maxPeak+30;

% plot waveforms
subplot(2,2,1), plot(wvf(1,:),wvf(2,:));
hold on
xlim([0 1023]);
ylim([-PlotmaxPeak PlotmaxPeak]);
title('wire 1', 'fontSize',14,'fontweight','bold');
subplot(2,2,2), plot(wvf(1,:),wvf(3,:));
xlim([0 1023]);
ylim([-PlotmaxPeak PlotmaxPeak]);
title('wire 2', 'fontSize',14,'fontweight','bold');
subplot(2,2,3), plot(wvf(1,:),wvf(4,:));
xlim([0 1023]);
ylim([-PlotmaxPeak PlotmaxPeak]);
title('wire 3', 'fontSize',14,'fontweight','bold');
subplot(2,2,4), plot(wvf(1,:),wvf(5,:));
xlim([0 1023]);
ylim([-PlotmaxPeak PlotmaxPeak]);
title('wire 4', 'fontSize',14,'fontweight','bold');

keep('Dur', 'meanDur', 'maxPeak');







