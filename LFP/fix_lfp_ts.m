% fix_lfp_ts
% fix lfp timestamps for sessions that have recording starts and stops 

% check sessions in the hippocampus
% df = readtable('D:/ryanh/github/harvey_et_al_2020/Rdata_pae_track_cylinder_all_cells.csv');
% 
% sessions_from_df = unique(df.session);

data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData';
sessions = dir(fullfile(data_path,'*.mat'));
sessions = struct2table(sessions);
sessions = sessions.name;

for i = 1:length(sessions)
    data = load(fullfile(data_path,sessions{i}),'lfp');
    
    generated_ts = (0:length(data.lfp.ts)) / 200;
    
    last_frame_gen_ts(i) = generated_ts(end);
    last_lfp_ts(i) = data.lfp.ts(end);
    disp(sessions(i).name)
end

offset = abs(last_frame_gen_ts - last_lfp_ts);
list_of_violations = find(offset>1);


%% create lfp timestamps to load in cases where lfp timestamps are misaligned
WaitMessage = parfor_wait(length(list_of_violations),'Waitbar',false,'ReportInterval',1);

for i = list_of_violations
    load(fullfile(data_path,sessions{i}),'session_path','basename')
    csc_list = dir(fullfile(session_path,'*.ncs'));
    
    fn = fullfile(session_path,csc_list(1).name);
    ts = Nlx2MatCSC(fn,[1 0 0 0 0], 0, 1, [] );
    ts = (ts-ts(1)) ./ 10^6;

    lfp = bz_GetLFP(1,'basepath',session_path,...
        'basename',basename,...
        'noPrompts',true,...
        'downsample',1);

    corrected_ts = interp1(linspace(1,length(lfp.data),length(ts)), ts, 1:length(lfp.data));
    
    writeNPY(corrected_ts, fullfile(session_path,'lfp_ts.npy'))
    
    WaitMessage.Send;
end
WaitMessage.Destroy

%%

