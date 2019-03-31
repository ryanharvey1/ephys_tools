function  [Start_ts, End_ts] =  FindRippleIntervals(eeg_tsd, threshold)
%   
%      [Start_ts, End_ts] =  FindRippleIntervals(eeg_tsd, threshold)
%
%   Find time intervals in which the EEG shows 'Ripples', i.e. intervals
%   in which the 100-300Hz filtered EEG signal amplitude is above 'threshold'.
%
%  Uses the envelope of absolute magnitude of the ripple-filtered eeg amplitude 
%  and determines the Start,Peak and End times of intervals in which the abs(filtered_amplitude)
%  is above the input argument 'threshold'.
%
%  IN: eeg_tsd       ...  tsd of FILTERED (100-300Hz) eeg
%      threshold     ...  Ripple interval begins when EEG enveloppe goes above threshold and
%                         ends when enveloppe sinks below threshold; 
%
% OUT: Start_ts, End_ts ... ts objects containing the start and end timestamps of a series 
%                           of intervals. Can be used directly as arguments for 
%                           the tsd/Restrict method
%
% PL  jan 2000
%

MinGapLength = 1000;   % minimum gap length  between valid ripples (0.1 sec)
MinRippleInterval = 10; % minimum length of a ripple interval      (1 msec)

abs_eeg = abs(Data(eeg_tsd));
ts_eeg  = Range(eeg_tsd,'ts');
[ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts_eeg,abs_eeg);

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

% make figure for visual inspection
figure;
%plot(ts_eeg,abs_eeg);
%axis tight;
%hold on;
plot(ts_env,pow_env,'r');
axis tight;
hold on;
startend_ts = [Data(Start_ts) Data(End_ts)];
plot(startend_ts',threshold*ones(size(startend_ts))');
hold off;
