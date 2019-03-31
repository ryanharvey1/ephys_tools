function  [Start_ts, End_ts] =  RippleIntervals(eeg_tsd, threshold, figure_flag, epochs, xlims)

% RippleIntervals  Finds (& plots if desired) time intervals in which EEG shows ripples
%   
% [Start_ts, End_ts] =  RippleIntervals(eeg_tsd, threshold, figure_flag, epochs, xlims)
%
% INPUTS: 
%   eeg_tsd - tsd of FILTERED (100-300Hz) eeg
%   threshold - Ripple interval begins when EEG envelope goes above threshold and
%       ends when envelope sinks below threshold; 
%   figure_flag - 1=plot a figure with visible threshold and intervals to check
%       visually the threshold setting, 0=don't plot  
%   epochs - struct with names and time intervals of epochs
%   xlims - x-axis range for graph - 1x2 array [start_timestamp end_timestamp]
% OUTPUTS: 
%   Start_ts, End_ts - ts objects containing the start and end timestamps of a series 
%       of intervals. Can be used directly as arguments for the tsd/Restrict method
%
% Find time intervals in which the EEG shows 'Ripples', i.e. intervals
%   in which the 100-300Hz filtered EEG signal amplitude is above 'threshold'.
% Uses the envelope of absolute magnitude of the ripple-filtered eeg amplitude 
%   and determines the Start,Peak and End times of intervals in which the 
%   abs(filtered_amplitude) is above the input argument 'threshold'.
%
% PL 1/00, last modified '04 by MN


if (nargin == 2)
    figure_flag = 0;
end


MinGapLength = 1000;   % minimum gap length  between valid ripples (0.1 sec)
MinRippleInterval = 10; % minimum length of a ripple interval      (1 msec)

abs_eeg = abs(Data(eeg_tsd));
ts_eeg  = Range(eeg_tsd,'ts');
clear eeg_tsd;
[ts_env, pow_env, PeaksIdx] = PowerEnvelope(ts_eeg,abs_eeg);
clear ts_eeg abs_eeg;

% find intervals with enveloppe above threshold
inth = find(pow_env > threshold);
ii   = find(diff(inth) > 2);
start_i = [inth(1); inth(ii+1)];
end_i = [inth(ii); inth(end)];

S = ts_env(start_i);
E = ts_env(end_i);

% cut out gaps shorter than MinGapLength
i = 2;
while i < length(S)
  while (S(i) - E(i-1)) < MinGapLength & i < length(S)
    S(i) = [];
    E(i-1) = [];
  end
  i = i+1;
end

% cut out ripple intervals shorter than MinRippleInterval
i=1;
while i <= length(S)
  while (E(i)-S(i)) < MinRippleInterval & i < length(S)
    S(i) = [];
    E(i) = [];
  end
  i = i+1;
end

Start_ts = ts(S);
End_ts = ts(E);

if (figure_flag > 0) 
    % make figure for visual inspection
    plot(ts_env/10000/60, pow_env, 'r');
    axis tight;
    hold on;
    startend_ts = [Data(Start_ts) Data(End_ts)];
    plot(startend_ts'/10000/60,threshold*ones(size(startend_ts))');
    text(-.12,.1,'Power envelope of','Units','normalized','FontSize',9,'Rotation',90);
    text(-.09,-.05,'abs(ripple-filtered eeg amplitude)','Units','normalized','FontSize',8,'Rotation',90);
    set(gca, 'XTick', []);

    
    t0 = xlims(1)/600000;  
    t1 = xlims(2)/600000;
    maxpowenv = max(pow_env);
    maxinterval = max(threshold*ones(size(startend_ts))');
    y1 = max([maxpowenv(1) maxinterval(1)]);
    y0 = 0;    
    DrawEpochLines(y1,y0,epochs);
    axis([t0 t1 y0 y1]);
    
end%if


%--------------------------------------------
function DrawEpochLines(y0,y1,epochs)

    for iep = 1:length(epochs.names)
        t0 = epochs.intervals{iep}(1)/10000/60;
        t1 = epochs.intervals{iep}(2)/10000/60;
        tm = (t0+t1)/2;
        hhl = line([t0 t0], [y0 y1]);
        set(hhl,'Color','r');
        hhl = line([t1 t1], [y0 y1]);
        set(hhl,'Color','b');
        text(tm,y0-(y1-y0)/15,epochs.names{iep});
    end