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

if isstruct(data)
    % unpack input
    p = inputParser;
    p.addParameter('figs',0);
    p.addParameter('lfp',data.lfp.signal);
    p.addParameter('lfp_ts',data.lfp.ts);
    p.addParameter('frequency',data.lfp.lfpsamplerate);
    p.addParameter('speed',data.frames(:,5));
    p.addParameter('mov_ts',data.frames(:,1));
    p.addParameter('ripple_fs',[130 200]);
    p.addParameter('save',1);
    p.addParameter('overwrite',0);

    p.parse(varargin{:});
    figs = p.Results.figs;
    lfp = p.Results.lfp;
    lfp_ts = p.Results.lfp_ts;
    frequency = p.Results.frequency;
    speed = p.Results.speed;
    mov_ts = p.Results.mov_ts;
    passband = p.Results.ripple_fs;
    save_results = p.Results.save;
    overwrite = p.Results.overwrite;

elseif isempty(data)
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
    
    if all_sessions
        data = load(fullfile(sessions(s).folder,sessions(s).name));
        
        p = inputParser;
        p.addParameter('figs',0);
        p.addParameter('lfp',data.lfp.signal);
        p.addParameter('lfp_ts',data.lfp.ts);
        p.addParameter('frequency',data.lfp.lfpsamplerate);
        p.addParameter('speed',data.frames(:,5));
        p.addParameter('mov_ts',data.frames(:,1));
        p.addParameter('ripple_fs',[130 200]);
        p.addParameter('save',1);
        p.addParameter('overwrite',0);
        
        p.parse(varargin{:});
        figs = p.Results.figs;
        lfp = p.Results.lfp;
        lfp_ts = p.Results.lfp_ts;
        frequency = p.Results.frequency;
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
    for ch = 1:size(lfp,1)
        signal_filtered(ch,:) = BandpassFilter(lfp(ch,:), frequency, passband);
        
        r_ripple(ch) = snr(lfp_ts,signal_filtered(ch,:));

        r(ch) = snr(lfp_ts,lfp(ch,:));
    end
    r(isinf(r)) = NaN;
        
    % check if max SNR is < 1
    if nanmax(r) < 1
        ripple_info = NaN;
        disp('Recording was too noisy')
        if all_sessions
            clearvars -except varargin all_sessions sessions s sess_list ...
                save overwrite WaitMessage
            continue
        else
            return
        end
    end
        
    [~,good_channel] = nanmax(r);
    if nanmin(r) < 1
        [~,noise_channel] = nanmin(r);
        noise = lfp(noise_channel,:)';
    else
        noise_channel = NaN;
        noise = zeros(size(lfp,2),1);
    end
    
    signal = lfp(good_channel,:);
    
    % detect ripples
    [ripples] = bz_FindRipples_ephys_tools(signal',lfp_ts',...
        'EMGfilename',emg_file,...
        'EMGThresh',0.9,...
        'noise',noise,...
        'passband',passband,...
        'frequency',frequency);
    
    % store channel used, noise channel, and snr
    ripples.detectorinfo.detectionparms.channel_used = good_channel;
    ripples.detectorinfo.detectionparms.noise_channel = noise_channel;
    ripples.detectorinfo.detectionparms.snr = r;
    ripples.detectorinfo.detectionparms.snr_ripple_fs = r_ripple;

    % exclude movement
    temp_speed = interp1(mov_ts,speed,ripples.timestamps);
    above_speed_thres = ~any(temp_speed < 5 | isnan(temp_speed),2);
    
    ripples.timestamps(above_speed_thres,:) = [];
    ripples.peaks(above_speed_thres) = [];
    ripples.peakNormedPower(above_speed_thres,:) = [];
    
    disp(['After speed thresholding: ' num2str(length(ripples.peaks)) ' events.']);
    
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
    
    [maps,ripple_data,stats] = bz_RippleStats(signal_filtered(good_channel,:)',...
        lfp_ts',ripples,'frequency',frequency);
    
    [unfiltered,~,~] = bz_RippleStats(signal',lfp_ts',ripples,'frequency',frequency);
    maps.unfiltered_ripples = unfiltered.ripples;
    
    
    % pull out unfiltered and filtered ripples from lfp
    for r = 1:size(ripples.timestamps,1)
        idx = lfp_ts >= ripples.timestamps(r,1) & lfp_ts <= ripples.timestamps(r,2);
        ripples.unfiltered_ripple{r} = signal(idx);
        ripples.filtered_ripple{r} = signal_filtered(good_channel,idx);
    end
    
    ripple_info.ripples = ripples;
    ripple_info.maps = maps;
    ripple_info.ripple_data = ripple_data;
    ripple_info.stats = stats;
    ripple_info.ripples.detectorinfo.ProcessedDatafile =...
        fullfile(sessions(s).folder,sessions(s).name);
    
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
