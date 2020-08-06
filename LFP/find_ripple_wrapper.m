function ripple_info = find_ripple_wrapper(data,varargin)
%find_ripple_wrapper: wrapper for findRipples from TSToolbox_Utils
%
%   Input:
%           data: ephys_tools data structure. (if you don't have this, you
%                                               can use [] for the first input,
%                                               also this can be your projects
%                                               base dir to run through all sessions)
%           varargin:
%               figs: debugging figs, 0 or 1
%               lfp: channels x time
%               lfp_ts: time stamps from lfp in seconds
%               speed: speed of animal (cm/sec)
%               mov_ts: time stamps from movement in seconds
%               save: 0 or 1, default = 1, if using ephys_tools
%               overwrite: 0 or 1, default = 0, if using ephys_tools
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
% dependent on buzcode.  ephys_tools\external_packages\buzcode
%
% Ryan Harvey 2019

sess_list = 1;
all_sessions = 0;

if isstruct(data) % if the user supplied the ephys_tools data structure

    lfp = bz_GetLFP('all','basepath',data.session_path,...
        'basename',data.basename,...
        'noPrompts',true,...
        'downsample',1);
    
    % unpack input
    p = inputParser;
    p.addParameter('figs',0);
    p.addParameter('speed',data.frames(:,5));
    p.addParameter('mov_ts',data.frames(:,1));
    p.addParameter('ripple_fs',[130 200]);
    p.addParameter('save',1);
    p.addParameter('overwrite',0);

    p.parse(varargin{:});
    figs = p.Results.figs;
    lfp_ts = lfp.timestamps';
    frequency = lfp.samplingRate;
    lfp = double(lfp.data');
    speed = p.Results.speed;
    mov_ts = p.Results.mov_ts;
    passband = p.Results.ripple_fs;
    save_results = p.Results.save;
    overwrite = p.Results.overwrite;

elseif isempty(data) % if the user supplied each variable individually
    p = inputParser;
    p.addParameter('figs',0);
    p.addParameter('lfp',[]);
    p.addParameter('lfp_ts',[]);
    p.addParameter('frequency',[]);
    p.addParameter('speed',[]);
    p.addParameter('mov_ts',[]);
    p.addParameter('ripple_fs',[130 200]);
    
    p.parse(varargin{:});
    figs = p.Results.figs;
    lfp = p.Results.lfp;
    lfp_ts = p.Results.lfp_ts;
    frequency = p.Results.frequency;
    speed = p.Results.speed;
    mov_ts = p.Results.mov_ts;
    passband = p.Results.ripple_fs;
    
elseif ischar(data) % can loop through all sessions in basedir
    p = inputParser;
    p.addParameter('save',1);
    p.addParameter('overwrite',0);
    p.parse(varargin{:});
    save_results = p.Results.save;
    overwrite = p.Results.overwrite;
    
    all_sessions = 1;
    
    basepath = data;
    
    % set up folder for data
    if ~exist(fullfile(basepath,'SWR'),'dir')
        mkdir(fullfile(basepath,'SWR'))
    end
    
    % get all sessions
    sessions = dir(fullfile(basepath,'ProcessedData','*.mat'));
    
    % check sessions that have already run
    finished_sess = dir(fullfile(basepath,'SWR','*.mat'));
    
    sess_list = find(~contains(extractBefore({sessions.name},'.mat'),...
        extractBefore({finished_sess.name},'_swr')));
    
    if overwrite
        sess_list = 1:length(sessions);
    end
end

if all_sessions
    WaitMessage = parfor_wait(length(sess_list),'Waitbar',true);
end

for s = sess_list
    if all_sessions % to loop through a data set
        data = load(fullfile(sessions(s).folder,sessions(s).name),...
            'frames','session_path','basename','sessionID','rat',...
            'mazetypes','linear_track','events','lfp','Spikes');
        lfp = bz_GetLFP('all','basepath',data.session_path,...
            'basename',data.basename,...
            'noPrompts',true,...
            'downsample',1);

        p = inputParser;
        p.addParameter('figs',0);
        p.addParameter('speed',data.frames(:,5));
        p.addParameter('mov_ts',data.frames(:,1));
        p.addParameter('ripple_fs',[130 200]);
        p.addParameter('save',1);
        p.addParameter('overwrite',0);
        
        p.parse(varargin{:});
        figs = p.Results.figs;
        lfp_ts = lfp.timestamps';
        frequency = lfp.samplingRate;
        lfp = double(lfp.data');
        speed = p.Results.speed;
        mov_ts = p.Results.mov_ts;
        passband = p.Results.ripple_fs;
        save_results = p.Results.save;
        overwrite = p.Results.overwrite;
    end
    
    processedpath=strsplit(data.session_path,filesep);
    processedpath(end-2:end)=[];
    emg_file = fullfile(strjoin(processedpath,filesep),'EMG_from_LFP',...
        [data.rat,'_',data.sessionID,'_emg.mat']);
            
    % find good and bad channels based on signal to noise ratio
    [signal_filtered,r_ripple,r] = get_signal_to_noise(lfp,lfp_ts,frequency,passband);
    
    [good_channel,signal,noise,noise_channel] = find_good_bad_channel(data,lfp,r);

    for ch = 1:size(lfp,1)
        if ~any(lfp(ch,:))
            continue
        end
        if length(data.lfp.channel_list.tetrode_num) > 8 % only use EMG on hyperdrives
            [ripples{ch}] = bz_FindRipples_ephys_tools(lfp(ch,:)',lfp_ts',...
                'EMGfilename',emg_file,...
                'EMGThresh',0.9,...
                'noise',noise,...
                'passband',passband,...
                'frequency',frequency);
        else
            [ripples{ch}] = bz_FindRipples_ephys_tools(lfp(ch,:)',lfp_ts',...
                'noise',noise,...
                'passband',passband,...
                'frequency',frequency);
        end
    end
    
    [ripples] = combine_and_exclude_close_events(ripples);
    
%     ripples = exclude_by_unit_activity(ripples,data,lfp_ts);
%     ripples = exclude_by_multi_unit_activity(ripples,data,lfp_ts);
    
    % store channel used, noise channel, and snr
    ripples.detectorinfo.detectionparms.channel_used = good_channel;
    ripples.detectorinfo.detectionparms.noise_channel = noise_channel;
    ripples.detectorinfo.detectionparms.snr = r;
    ripples.detectorinfo.detectionparms.snr_ripple_fs = r_ripple;

    % exclude movement
    [mov_ts,idx] = unique(mov_ts);
    temp_speed = interp1(mov_ts,speed(idx),ripples.timestamps);
    above_speed_thres = ~any(temp_speed < 3 | isnan(temp_speed),2);
    ripples.timestamps(above_speed_thres,:) = [];
    ripples.peaks(above_speed_thres) = [];
    ripples.peakNormedPower(above_speed_thres,:) = [];
    ripples.ch_map(above_speed_thres,:) = [];
    disp(['After speed thresholding: ' num2str(length(ripples.peaks)) ' events.']);
    
    % remove channels with < 2 ripples
    ripples = remove_ch_with_few_ripples(ripples);
    
    % if less than 2 ripples are found
    if size(ripples.timestamps,1) < 2
        ripple_info = NaN;
        disp('No ripples')
        if all_sessions
            clearvars -except varargin all_sessions sessions s sess_list ...
                save overwrite WaitMessage
            continue
        else
            return
        end
    end
    
    [maps,ripple_data,stats] = get_ripple_stats(lfp,signal_filtered,...
        ripples,frequency,lfp_ts);
    
    % pull out unfiltered and filtered ripples from lfp
    for r = 1:size(ripples.timestamps,1)
        idx = lfp_ts >= ripples.timestamps(r,1) & lfp_ts <= ripples.timestamps(r,2);
        ripples.unfiltered_ripple{r} = lfp(ripples.ch_map(r),idx);
        ripples.filtered_ripple{r} = signal_filtered(ripples.ch_map(r),idx);
    end
    
    ripple_info.ripples = ripples;
    ripple_info.maps = maps;
    ripple_info.ripple_data = ripple_data;
    ripple_info.stats = stats;
    if exist('sessions','var')
        ripple_info.ripples.detectorinfo.ProcessedDatafile =...
            fullfile(sessions(s).folder,sessions(s).name);
    end
    
    if figs
        bz_PlotRippleStats(ripple_info.maps,ripple_info.ripple_data,...
            ripple_info.stats,'frequency',frequency);
        ripple_figs.lfp_viewer(ripple_info,data)
        ripple_figs.ripple_grid(ripple_info)
        ripple_figs.ripple_per_channel(ripple_info,data)
        ripple_figs.ripple_location(ripple_info,data)
    end

    if all_sessions || save_results && overwrite
        % save mat file to processed data folder
        processedpath=strsplit(data.session_path,filesep);
        processedpath(end-2:end)=[];
        save(fullfile(strjoin(processedpath,filesep),'SWR',...
            [data.rat,'_',data.sessionID,'_swr']),'-struct','ripple_info','-v7.3')
        
        clearvars -except varargin all_sessions sessions s sess_list save ...
            overwrite WaitMessage
        WaitMessage.Send;
    end
end
if all_sessions
    WaitMessage.Destroy;
end
end

function [maps,ripple_data,stats] = get_ripple_stats(lfp,signal_filtered,ripples,frequency,lfp_ts)
ripples_ = [];
frequency_ = [];
phase_ = [];
amplitude = [];
unfiltered_ripples = [];
peakFrequency = [];
peakAmplitude = [];
duration_ = [];

for ch = 1:size(lfp,1)
    idx = ripples.ch_map == ch;
    if ~any(idx)
        continue
    end
    temp_ripples = ripples;
    temp_ripples.peaks(~idx) = [];
    temp_ripples.timestamps(~idx,:) = [];
    temp_ripples.peakNormedPower(~idx) = [];
    
    [maps_,ripple_data_,~] = bz_RippleStats(signal_filtered(ch,:)',...
        lfp_ts',temp_ripples,'frequency',frequency);
    
    [unfiltered,~,~] = bz_RippleStats(lfp(ch,:)',lfp_ts',temp_ripples,'frequency',frequency);
    maps_.unfiltered_ripples = unfiltered.ripples;
    
    if size(maps_.ripples,2) ~= size(ripples_,2) && ~isempty(ripples_)
        maps_ = sync_to_size(maps_,size(ripples_,2));
    end
        
    ripples_ = [ripples_;maps_.ripples];
    frequency_ = [frequency_;maps_.frequency];
    phase_ = [phase_;maps_.phase];
    amplitude = [amplitude;maps_.amplitude];
    unfiltered_ripples = [unfiltered_ripples;maps_.unfiltered_ripples];
    
    peakFrequency = [peakFrequency;ripple_data_.peakFrequency];
    peakAmplitude = [peakAmplitude;ripple_data_.peakAmplitude];
    duration_ = [duration_;ripple_data_.duration];
    
end

corrBinSize = 0.01;
[stats.acg.data,stats.acg.t] = CCG(ripples.peaks,ones(length(ripples.peaks),1),'binSize',corrBinSize);
[stats.amplitudeFrequency.rho,stats.amplitudeFrequency.p] = corrcoef(peakAmplitude,peakFrequency);
[stats.durationFrequency.rho,stats.durationFrequency.p] = corrcoef(duration_,peakFrequency);
[stats.durationAmplitude.rho,stats.durationAmplitude.p] = corrcoef(duration_,peakAmplitude);

maps.ripples = ripples_;
maps.frequency = frequency_;
maps.phase = phase_;
maps.amplitude = amplitude;
maps.unfiltered_ripples = unfiltered_ripples;

ripple_data.peakFrequency = peakFrequency;
ripple_data.peakAmplitude = peakAmplitude;
ripple_data.duration = duration_;

end

function maps = sync_to_size(maps,bins)
for i = 1:size(maps.ripples,1)
    
    ripples(i,:) = interp1(linspace(1,bins,size(maps.ripples,2)),...
        maps.ripples(i,:),1:bins);
    
    frequency(i,:) = interp1(linspace(1,bins,size(maps.frequency,2)),...
        maps.frequency(i,:),1:bins);
    
    phase(i,:) = interp1(linspace(1,bins,size(maps.phase,2)),...
        maps.phase(i,:),1:bins);
    
    amplitude(i,:) = interp1(linspace(1,bins,size(maps.amplitude,2)),...
        maps.amplitude(i,:),1:bins);
    
    unfiltered_ripples(i,:) = interp1(linspace(1,bins,size(maps.unfiltered_ripples,2)),...
        maps.unfiltered_ripples(i,:),1:bins);
end
maps.ripples = ripples;
maps.frequency = frequency;
maps.phase = phase;
maps.amplitude = amplitude;
maps.unfiltered_ripples = unfiltered_ripples;
end

function ripples = exclude_by_unit_activity(ripples,data,lfp_ts)
dt = (lfp_ts(2) - lfp_ts(1));
bin_edges = [lfp_ts - dt/2, lfp_ts(end) + dt/2];
for i = 1:length(data.Spikes)
    binned_spikes(i,:) = histcounts(data.Spikes{i},bin_edges);
end
% conver to binary matrix
binned_spikes = binned_spikes > 0;

for r = 1:size(ripples.timestamps,1)
    idx = lfp_ts >= ripples.timestamps(r,1) & lfp_ts <= ripples.timestamps(r,2);
    sum_spikes(:,r) = sum(binned_spikes(:,idx),2);
end
% convert to binary
sum_spikes = sum_spikes > 0;
idx = sum(sum_spikes) < 1;

ripples.peaks(logical(idx)) = [];
ripples.timestamps(logical(idx),:) = [];
ripples.peakNormedPower(logical(idx)) = [];
ripples.ch_map(logical(idx)) = [];

% idx = lfp_ts >= ripples.timestamps(116,1) & lfp_ts <= ripples.timestamps(116,2);
% 
% data.lfp.signal(:,idx)
end

function ripples = exclude_by_multi_unit_activity(ripples,data,lfp_ts)
minimum_duration = 0.015;
zscore_thres = 0;
bin_size = 0.05;

ts_ifr = 0:bin_size:lfp_ts(end);
for k = 1:length(data.Spikes)
    [ifr(:,k),~]=instantfr(data.Spikes{k},ts_ifr);
end

% zscore and mean
zscored_data = mean(zscore(ifr),2);
% threshold
is_above_mean = zscored_data >= zscore_thres;
% find epochs above threshold
indexout=contiguousframes(is_above_mean,floor(bin_size/minimum_duration));
[start,ends,ngroups]=findgroups(indexout);

% get all timestamps in range
times = [];
for i = 1:ngroups
    times = [times;ts_ifr(start(i):ends(i))'];
end

% compare ripple peak time against mua activity
for r = 1:size(ripples.timestamps,1)
    if min(abs(ripples.peaks(r) - times)) <= minimum_duration * 2
        idx(r) = 1;
    else
        idx(r) = 0;
    end
end

ripples.peaks(logical(idx)) = [];
ripples.timestamps(logical(idx),:) = [];
ripples.peakNormedPower(logical(idx)) = [];
ripples.ch_map(logical(idx)) = [];
end

function ripples = combine_and_exclude_close_events(ripples)
% unpack ripple events on each channel & exclude successive events
% that occur within a `close_event_threshold` of a previously occuring event.

% pull out data across channels
peak_ = [];
ch_map = [];
timestamps = [];
peakNormedPower = [];
for ch = 1:size(ripples,2)
    if isempty(ripples{ch})
        continue
    end
    peak_ = [peak_;ripples{ch}.peaks];
    timestamps = [timestamps;ripples{ch}.timestamps];
    peakNormedPower = [peakNormedPower;ripples{ch}.peakNormedPower];
    ch_map = [ch_map;repmat(ch,length(ripples{ch}.peaks),1)];
    stdev(ch) = ripples{ch}.stdev;
end
detectionparms = ripples{ch}.detectorinfo.detectionparms;
detectionparms.lfp = [];

% sort
[peak_,idx] = sort(peak_);
candidate_event_times = timestamps(idx,:);
ch_map = ch_map(idx,:);
peakNormedPower = peakNormedPower(idx,:);

% find close ripples across channels and remove them
close_event_threshold = 0.00;
n_events = size(candidate_event_times,1);
new_event_index = [1:n_events]';
new_event_times = candidate_event_times;
for ind = 1:n_events
    is_too_close = candidate_event_times(ind,2) + close_event_threshold >...
        new_event_times(:,1) & new_event_index > ind;
    new_event_index = new_event_index(~is_too_close);
    new_event_times = new_event_times(~is_too_close,:);
end
peak_ = peak_(new_event_index);
peakNormedPower = peakNormedPower(new_event_index);
ch_map = ch_map(new_event_index);

% pull out some metadata
detectorname = ripples{1}.detectorinfo.detectorname;
detectiondate = ripples{1}.detectorinfo.detectiondate;
detectionintervals = ripples{1}.detectorinfo.detectionintervals;

% recreate ripples struct
ripples = [];
ripples.peaks = peak_;
ripples.timestamps = new_event_times;
ripples.peakNormedPower = peakNormedPower;
ripples.ch_map = ch_map;
ripples.stdev = stdev;
ripples.detectorinfo.detectionparms = detectionparms;
ripples.detectorinfo.detectorname = detectorname;
ripples.detectorinfo.detectiondate = detectiondate;
ripples.detectorinfo.detectionintervals = detectionintervals;

% remove channels with < 2 ripples
ripples = remove_ch_with_few_ripples(ripples);
end

function ripples = remove_ch_with_few_ripples(ripples)
% remove channels with < 2 events
idx = zeros(length(ripples.ch_map),1);
for ch = unique(ripples.ch_map)'
    if sum(ripples.ch_map == ch) < 2
        idx(ripples.ch_map == ch) = 1;
    end
end
ripples.peaks(logical(idx)) = [];
ripples.timestamps(logical(idx),:) = [];
ripples.peakNormedPower(logical(idx)) = [];
ripples.ch_map(logical(idx)) = [];
end

function [signal_filtered,r_ripple,r] = get_signal_to_noise(lfp,lfp_ts,frequency,passband)

for ch = 1:size(lfp,1)
    signal_filtered(ch,:) = BandpassFilter(lfp(ch,:), frequency, passband);
    
    r_ripple(ch) = snr(lfp_ts,signal_filtered(ch,:));
    
    r(ch) = snr(lfp_ts,lfp(ch,:));
end
r(isinf(r)) = NaN;
r_ripple(isinf(r_ripple)) = NaN;
end

function [good_channel,signal,noise,noise_channel] = find_good_bad_channel(data,lfp,r)

% find disconnected channels
session_info = LoadParameters(data.session_path);
% mark good channels
good = zeros(1,size(lfp,1));
% check spike group from xml file
good([session_info.spikeGroups.groups{:}] + 1) = 1;
% check if a channel is all zeros
connected_channels(sum(lfp') ~= 0 & good,1) = 1;

% pick channel to use
good = find(connected_channels);
[~,I] = nanmax(r(good));
good_channel = good(I);

signal = lfp(good_channel,:);

% pick noise channel based on xml 
if any(connected_channels == 0)
    noise_channel = find(connected_channels == 0,1,'first');
    noise = lfp(noise_channel,:)';
% elseif nanmin(r) < 1 && nanmax(r) > 1
%     [~,noise_channel] = nanmin(r);
%     noise = lfp(noise_channel,:)';
else
    noise_channel = NaN;
    noise = zeros(size(lfp,2),1);
end
end