function ripple_info = find_ripple_wrapper(data,varargin)
%find_ripple_wrapper: wrapper for findRipples from TSToolbox_Utils
%
%   Input:
%           data: ephys_tools data structure. (if you don't have this, you 
%                                               can use [] for the first input)
%           varargin:
%               figs: debugging figs, 0 or 1
%               lfp: channels x time 
%               lfp_ts: time stamps from lfp in seconds
%               speed: speed of animal (cm/sec)
%               mov_ts: time stamps from movement in seconds
%
%   Output:
%           ripple_info:
%                ripple_info.ripple_start: start time of ripple
%                ripple_info.ripple_peak: peak time of ripple
%                ripple_info.ripple_end: end time of ripple
%                ripple_info.peakNormalizedPower: power of ripple
%                ripple_info.frequency: frequency of ripple
%                ripple_info.amplitude: peak amplitude of ripple
%                ripple_info.ripple_dur: duration of ripple in seconds
%                ripple_info.unfiltered_ripple: unfiltered ripple
%                ripple_info.filtered_ripple: filtered ripple
%                ripple_info.time: time stamps during ripple
%
% Ryan Harvey 2019

if isstruct(data)
    % unpack input
    p = inputParser;
    p.addParameter('figs',0);
    p.addParameter('lfp',data.lfp.signal);
    p.addParameter('lfp_ts',data.lfp.ts);
    p.addParameter('frequency',data.lfp.lfpsamplerate);
    p.addParameter('speed',data.frames(:,5));
    p.addParameter('mov_ts',data.frames(:,1));
    
    p.parse(varargin{:});
    figs = p.Results.figs;
    lfp = p.Results.lfp;
    lfp_ts = p.Results.lfp_ts;
    frequency = p.Results.frequency;
    speed = p.Results.speed;
    mov_ts = p.Results.mov_ts;
else
    p = inputParser;
    p.addParameter('figs',0);
    p.addParameter('lfp',[]);
    p.addParameter('lfp_ts',[]);
    p.addParameter('frequency',[]);
    p.addParameter('speed',[]);
    p.addParameter('mov_ts',[]);
    
    p.parse(varargin{:});
    figs = p.Results.figs;
    lfp = p.Results.lfp;
    lfp_ts = p.Results.lfp_ts;
    frequency = p.Results.frequency;
    speed = p.Results.speed;
    mov_ts = p.Results.mov_ts;
end

for ch =1:size(lfp,1)
    signal = lfp(ch,:);
    
    % filter signal to ripple range [150 250 hz]
    signal_filtered = BandpassFilter(signal, 1000, [150 250]);
    
    % locate potential ripples
    [ripples,~] = findRipples(signal_filtered','frequency',frequency);
    
    if isempty(ripples)
        continue
    end
    
    % align to session
    ripples(:,1:3) = ripples(:,1:3) - lfp_ts(1);
    
    % exclude movement
    temp_speed = interp1(mov_ts,speed,ripples(:,1:3));
    ripples = ripples(any(temp_speed < 4 | isnan(temp_speed),2),:);
    
    % max amplitude
    for r = 1:size(ripples)
        ripples(r,6) = max(signal_filtered(lfp_ts >= ripples(r,1) & lfp_ts <= ripples(r,3)));
    end
    
    % number and mean power
    n_ripples(ch) = (size(ripples,1));
    mean_power(ch) = nanmean(ripples(:,4));
    
    % store ripple info
    all_ripples{ch} = ripples;
end

% if no ripple is found
if isempty(ripples)
    ripple_info.ripple_start = NaN;
    ripple_info.ripple_peak = NaN;
    ripple_info.ripple_end = NaN;
    ripple_info.peakNormalizedPower = NaN;
    ripple_info.frequency = NaN;
    ripple_info.amplitude = NaN;
    ripple_info.ripple_dur = NaN;
    ripple_info.unfiltered_ripple = NaN;
    ripple_info.filtered_ripple = NaN;
    ripple_info.time = NaN;
    return
end

% locate channel with highest mean power and number of ripples
[~,ch] = max(mean_power.*n_ripples);
ripples = all_ripples{ch};
signal = lfp(ch,:);

% re-bandpass filter
signal_filtered = BandpassFilter(signal, 1000, [150 250]);

% get ripple duration in sec
ripple_dur = ripples(:,3) - ripples(:,1);

% pull out unfiltered and filtered ripples from lfp
for r = 1:size(ripples,1)
    idx = lfp_ts >= ripples(r,1) & lfp_ts <= ripples(r,3);
    unfiltered_ripple{r} = signal(idx);
    filtered_ripple{r} = signal_filtered(idx);
    time{r} = lfp_ts(idx);
end

% pack into struct
ripple_info.ripple_start = ripples(:,1);
ripple_info.ripple_peak = ripples(:,2);
ripple_info.ripple_end = ripples(:,3);
ripple_info.peakNormalizedPower = ripples(:,4);
ripple_info.frequency = ripples(:,5);
ripple_info.amplitude = ripples(:,6);
ripple_info.ripple_dur = ripple_dur;
ripple_info.unfiltered_ripple = unfiltered_ripple;
ripple_info.filtered_ripple = filtered_ripple;
ripple_info.time = time;


if figs
    figure;
    p = ceil(sqrt(size(ripples,1)));
    colors = viridis(length(unique(ripples(:,4))));
    for r = 1:size(ripples,1)
        subplot(p,p,r)
        plot(time{r},unfiltered_ripple{r},'color',interp1(unique(ripples(:,4)),colors,ripples(r,4)))
        hold on
        box off
        axis off
    end
    sgtitle([num2str(r),' Ripples on Channel ',num2str(ch),', Color code: power'],'Color','w')
    darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
    
%     figure
%     for r = 1:size(ripples,1)
%         subplot(p,p,r)
%         idx = lfp_ts >= ripples(r,1)-0.5 & lfp_ts <= ripples(r,3)+0.5;
%         pspectrum(signal(idx),frequency,'spectrogram',...
%             'OverlapPercent',99,'MinTHreshold',10,'FrequencyLimits',[100, 250])
%         colorbar off
%         axis off
%         box off
%         title('')
%         colormap(viridis)
%         hold on
%         pause(.0001)
%     end
    
    figure;
    plot(lfp_ts,signal);
    hold on
    yLim = ylim;
    for j=1:size(ripples,1)
        plot([ripples(j,1) ripples(j,1)],yLim,'g-');
        plot([ripples(j,2) ripples(j,2)],yLim,'k-');
        plot([ripples(j,3) ripples(j,3)],yLim,'r-');
    end
end
end

