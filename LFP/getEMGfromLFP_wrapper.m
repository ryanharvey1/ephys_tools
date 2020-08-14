function getEMGfromLFP_wrapper(basepath)
% getEMGfromLFP_wrapper: run through all sessions and make emg file
% basepath is the path to your project
% Ryan H 2020

% set up folder for data
if ~exist(fullfile(basepath,'EMG_from_LFP'),'dir')
    mkdir(fullfile(basepath,'EMG_from_LFP'))
end

% get all sessions
sessions = dir(fullfile(basepath,'ProcessedData','*.mat'));

% check sessions that have already run
finished_sess = dir(fullfile(basepath,'EMG_from_LFP','*.mat'));

sess_list = find(~contains(extractBefore({sessions.name},'.mat'),...
    extractBefore({finished_sess.name},'_emg')));
sessions = sessions(sess_list);

function run_save_emg(data)
emg = getEMGfromLFP(data.lfp.signal,'fs',data.lfp.lfpsamplerate,'graphics',false);

% save mat file to processed data folder
processedpath=strsplit(data.session_path,filesep);
processedpath(end-2:end)=[];
save(fullfile(strjoin(processedpath,filesep),'EMG_from_LFP',...
    [data.rat,'_',data.sessionID,'_emg']),'-struct','emg','-v7.3')
end
fcn = @run_save_emg;


WaitMessage = parfor_wait(length(sess_list),'Waitbar',true);
parfor i = 1:length(sessions)
    
    data = load(fullfile(sessions(i).folder,sessions(i).name),...
        'session_path','rat','sessionID');
    
    [~,data.basename] = fileparts(data.session_path);
    
%     sess_info = LoadParameters(data.session_path);
    
    lfp = bz_GetLFP('all','basepath',data.session_path,...
        'basename',data.basename,...
        'noPrompts',true,...
        'downsample',1);
    data.lfp.signal = double(lfp.data);
    data.lfp.lfpsamplerate = lfp.samplingRate;

    feval(fcn, data); 
    
    WaitMessage.Send;
end
WaitMessage.Destroy;
end