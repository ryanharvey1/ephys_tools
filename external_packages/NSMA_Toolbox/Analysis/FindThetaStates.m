function  [Start_ts, End_ts] =  FindThetaStates(eeg_tsd, time_fraction)
%   
%      [Start_ts, End_ts] =  FindThetaStates(eeg_tsd, time_fraction)
%
%   Find time intervals in which the animal is presumably in 'theta state'.
%
%  Uses a crude local power estimate (envelope of theta-filtered eeg amplitude squared)
%  and determines a threshold so that the animal is a fraction time_fraction 
%  of the total time span in theta state. (the fraction of time is supposed to be determined
%  by fraction of time the animal spent 'in motion')
%
%  IN: eeg_tsd       ...  tsd of FILTERED (6-10Hz) eeg
%      time_fraction ... fraction of time (0 < time_fraction < 1) the animal is in theta; 
%                          this is equivalent to choosing a threshold for the local power
%                          enveloppe
%
% OUT: Start_ts, End_ts ... ts objects containing the start and end timestamps of a series 
%                           of intervals. Can be used directly as arguments for 
%                           the tsd/Restrict method
%
% PL  dec 99
%

MinGapLength = 10000;   % minimum gap length (LIA) between theta cycles (1 sec)
MinThetaInterval = 10000; % minimum length of a theta interval (1 sec)

pow = Data(eeg_tsd).^2;
ts_eeg  = Range(eeg_tsd,'ts');
[ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts_eeg,pow);

% calculate power threshold from given time_fraction
[pow_hist,xbins] = hist(pow_env,100);     
pow_cpdf = cumsum(pow_hist)/sum(pow_hist);
tail = find(pow_cpdf > 1-time_fraction);
pow_threshold = xbins(tail(1))

% find intervals with thetapower above threshold
inth = find(pow_env > pow_threshold);
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

% cut out theta intervals shorter than MinThetaInterval
i=1;
while i <= length(S)
  while (E(i)-S(i)) < MinThetaInterval & i < length(S)
    S(i) = [];
    E(i) = [];
  end
  i = i+1;
end

Start_ts = ts(S);
End_ts = ts(E);

% make figure for visual inspection
figure;
plot(ts_eeg,pow);
axis tight;
hold on;
plot(ts_env,pow_env,'r');
startend_ts = [Data(Start_ts) Data(End_ts)];
plot(startend_ts',pow_threshold*ones(size(startend_ts))');
hold off;
