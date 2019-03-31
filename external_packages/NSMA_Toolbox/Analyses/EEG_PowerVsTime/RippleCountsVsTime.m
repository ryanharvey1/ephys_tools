function RippleCountsVsTime(threshold, eeg_tsd, binsize_sec, epochs, xlims)

% RippleCountsVsTime  Plots ripple counts and ripple power over time
%
% RippleCountsVsTime(threshold, eeg_tsd, binsize_sec, epochs, xlims)
%
% INPUTS:
%   threshold - value of ripple-filtered eeg envelope above which a ripple is
%       considered to occur
%   eeg_tsd - tsd of eeg data from good theta channel
%   binsize_sec - size of timebins in which to compute avg. theta power
%   epochs - struct with names and time intervals of epochs
%   xlims - x-axis limits -- 1x2 array of start & end timestamps
% OUTPUTS:
%   (none)
%
% PL '02, last modified '04 by MN


% Chebychev bandpass filter in frequency range 100 300
Z = Filter4Ripples(eeg_tsd,100,300);   % Z is the NOT-SUBSAMPLED filtered eeg tsd (20-60Hz = Gamma, 100-300Hz = ripples)

tstart = StartTime(eeg_tsd);
tend = EndTime(eeg_tsd);

[Ripple_Start_ts, Ripple_End_ts] = RippleIntervals(Z, threshold, 0, epochs, xlims);  % threshold for ripple intervals needs to be adjusted visually; 100 is a good first try 
xbins = (StartTime(Z):binsize_sec*10000:EndTime(Z))';
Ripple_Mean_ts = (Data(Ripple_Start_ts) + Data(Ripple_End_ts))/2; 
Ripple_Counts = hist(Ripple_Mean_ts,xbins);

pow = Data(Z).^2;
ts_eeg  = Range(Z,'ts');
[ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts_eeg,pow);

% the power envelope is not equidistant in ts (ts is located at theta peaks!).
% interpolate with equidistant 100Hz sampling rate

pow_ts_interp = tstart:10:tend;            % 100 timestamp intervalls corresponds to 1000 Hz sampling rate
pow_env_interp = interp1(ts_env,pow_env,pow_ts_interp);

% binsize in number of sample points at 1000Hz sampling rate:
binsize = round(binsize_sec)*1000;

% truncate end so that length is exactly a  multiple of binsize
pow_ts_interp(end-mod(length(pow_ts_interp),binsize)+1:end) = [];  
pow_env_interp(end-mod(length(pow_env_interp),binsize)+1:end) = [];

%calculate averages over binsize by reshaping and averaging over columns
av_pow_ts = mean(reshape(pow_ts_interp',binsize,length(pow_ts_interp)/binsize))';
av_pow_env = mean(reshape(pow_env_interp',binsize,length(pow_env_interp)/binsize))';


% PLOT: av_pow_env vs. av_pow_ts and ripple counts
if ~isempty(av_pow_env)
    H = plot(xbins/10000/60,Ripple_Counts,'b', ...
        av_pow_ts/10000/60, av_pow_env*(max(Ripple_Counts)/max(av_pow_env)), 'g');
    
    if ~isempty(epochs)
        hold on
        t0 = xlims(1)/600000;
        t1 = xlims(2)/600000;
        y1 = max(Ripple_Counts);
        y0 = 0;    
        DrawEpochLines(y1,y0,epochs);
        axis([t0 t1 y0 y1]);
        hold off;
    end
    
    set(H(1), 'Color', [0 0.9 0]);
    set(H(2), 'LineWidth', 0.1);
    set(H(2), 'Color', [0.5 0.5 0.5]);
    xlabel('time (min)')
    %set(gca, 'YTick', []);
    str(1) = {'Ripple Counts and'};
    str(2) = {'avg. ripple power'};
    str(3) = {['in ' num2str(binsize_sec) ' sec bins']};
    h = ylabel(str,'Units','normalized');
    g=get(h,'Position');
    g(1) = g(1)-.045;
    ylabel(str,'Position',g,'Units','normalized');
end


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
