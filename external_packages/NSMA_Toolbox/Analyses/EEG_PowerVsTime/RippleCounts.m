function RippleCounts(eeg_tsd, ripple_start_ts, ripple_end_ts, epochs, xlims)

% RippleCounts  Plots counts of ripples in 5-sec bins vs. time
%
% RippleCounts(eeg_tsd, ripple_start_ts, ripple_end_ts, epochs, xlims)
%
% INPUTS:
%   eeg_tsd - tsd of FILTERED (100-300Hz) eeg
%   ripple_start_ts, ripple_end_ts - ts objects containing the start and end timestamps of 
%       a series of ripple intervals
%   epochs - struct with names and time intervals of epochs
%   xlims - x-axis range for graph - 1x2 array [start_timestamp end_timestamp]
% OUTPUTS:
%   (none)
%
% PL '02, last modified '04 by MN


xbins = (StartTime(eeg_tsd):5*10000:EndTime(eeg_tsd))';
binsize=5*1000;
pow = Data(eeg_tsd).^2;
ts_eeg  = Range(eeg_tsd,'ts');
tstart = StartTime(eeg_tsd);
tend = EndTime(eeg_tsd);
clear eeg_tsd; %pack;

[ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts_eeg,pow);
clear ts_eeg pow; %pack;

pow_ts_interp = tstart:10:tend;
pow_ts_interp(end-mod(length(pow_ts_interp),binsize)+1:end) = [];
av_pow_ts_ripple = mean(reshape(pow_ts_interp',binsize,length(pow_ts_interp)/binsize))';
%save temp_ripple;

pow_env_interp = interp1(ts_env,pow_env,pow_ts_interp); 
clear pow_ts_interp; %pack;

pow_env_interp(end-mod(length(pow_env_interp),binsize)+1:end) = [];
av_pow_env_ripple = mean(reshape(pow_env_interp',binsize,length(pow_env_interp)/binsize))';
%save temp_ripple av_pow_env_ripple av_pow_ts_ripple xbins;
%clear all; pack;

Ripple_Mean_ts = (Data(ripple_start_ts) + Data(ripple_end_ts))/2; 
Ripple_Counts = hist(Ripple_Mean_ts,xbins);

plot(xbins/10000/60,Ripple_Counts,'g');
xlabel('time (min)');
ylabel('Ripple Counts in 5-sec. bins');

hold on;
t0 = xlims(1)/600000;
t1 = xlims(2)/600000;
y1 = max(Ripple_Counts);
y0 = 0;    
DrawEpochLines(y1,y0,epochs);
axis([t0 t1 y0 y1]);
hold off;


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
        %text(tm,y0-(y1-y0)/15,epochs.names{iep});
    end
