%re_run_sessions
% scratch code I made to rerun sessions that were improperly post processed
% Ryan H
sessions=dir('D:\Projects\PAE_PlaceCell\ProcessedData\*.mat');
githubcommitdate='31-Mar-2019';

for i=1:size(sessions,1)
    load(fullfile(sessions(i).folder,sessions(i).name),'date_processed','session_path')
    if datenum(date_processed)>datenum(githubcommitdate)
        data=postprocess(session_path,120,'yes',0);
    end
end
