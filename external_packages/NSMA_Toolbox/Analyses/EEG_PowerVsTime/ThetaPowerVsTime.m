function ThetaPowerVsTime(eeg_tsd, binsize_sec, epochs, xlims)

% ThetaPowerVsTime  Plots theta power over time
%
% ThetaPowerVsTime(eeg_tsd, binsize_sec, epochs, xlims)
%
% INPUTS:
%   eeg_tsd - tsd of eeg data from good theta channel
%   binsize_sec - size of timebins in which to compute avg. theta power
%   epochs - struct with names and time intervals of epochs
%   xlims - x-axis limits -- 1x2 array of start & end timestamps
% OUTPUTS:
%   (none)
%
% PL '02, last modified '04 by MN


% filter for theta
[Y,Z] = Filter4Theta2(eeg_tsd,6,10);   % Z is the NOT-SUBSAMPLED filtered eeg tsd 

pow = Data(Z).^2;
ts_eeg  = Range(Z,'ts');
[ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts_eeg,pow);

% the power envelope is not equidistant in ts (ts is located at theta peaks!).
% interpolate with equidistant 100Hz sampling rate

tstart = StartTime(eeg_tsd);
tend = EndTime(eeg_tsd);

pow_ts_interp = tstart:100:tend;            % 100 timestamp intervals corresponds to 100 Hz sampling rate
pow_env_interp = interp1(ts_env,pow_env,pow_ts_interp);

% binsize in number of sample points at 100Hz sampling rate:
binsize = round(binsize_sec)*100;

% truncate end so that length is exactly a  multiple of binsize
pow_ts_interp(end-mod(length(pow_ts_interp),binsize)+1:end) = [];  
pow_env_interp(end-mod(length(pow_env_interp),binsize)+1:end) = [];

%calculate averages over binsize by reshaping and averaging over columns
av_pow_ts = mean(reshape(pow_ts_interp',binsize,length(pow_ts_interp)/binsize))';
av_pow_env = mean(reshape(pow_env_interp',binsize,length(pow_env_interp)/binsize))';


% PLOT: av_pow_env vs. av_pow_ts
H = plot(av_pow_ts/10000/60, av_pow_env);
hold on;
t0 = xlims(1)/600000;
t1 = xlims(2)/600000;
y1 = max(av_pow_env);
y0 = 0;    
% DrawEpochLines(y1,y0,epochs); Ryan H 12/9/16
axis([t0 t1 y0 y1]);
    
set(H(1), 'Color', [0.5 0.5 0.5]);
set(H(1), 'LineWidth', 0.1);
%xlabel('time (min)');
str(1) = {'Avg. \Theta power'};
str(2) = {['in ' num2str(binsize_sec) ' sec bins']};
h = ylabel(str,'Units','normalized');
g=get(h,'Position');
g(1) = g(1)-.025;
ylabel(str,'Position',g,'Units','normalized');


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
