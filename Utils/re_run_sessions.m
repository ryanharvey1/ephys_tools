%re_run_sessions
% scratch code I made to rerun sessions that were improperly post processed
% Ryan H
sessions=dir('D:\Projects\PAE_PlaceCell\ProcessedData\*.mat');
githubcommitdate='31-Mar-2019';

for i=1:size(sessions,1)
    load(fullfile(sessions(i).folder,sessions(i).name),'date_processed','session_path','mazetypes')
    if datenum(date_processed)>datenum(githubcommitdate) &&...
            any(contains(mazetypes,'cylinder','IgnoreCase',true)) ||...
            any(contains(mazetypes,'box','IgnoreCase',true))
        
        disp(session_path)
        
        track_length=TrackLength(session_path); % SET TRACK LENGTH
        
        data=postprocess(session_path,track_length,'yes',0);
    end
end


sessions=dir('D:\Projects\PAE_PlaceCell\ProcessedData\*.mat');
for i=1:size(sessions,1)
    load(fullfile(sessions(i).folder,sessions(i).name),'date_processed','session_path','mazetypes')
    if isempty(mazetypes{1})
        
        disp(session_path)
        
        track_length=TrackLength(session_path); % SET TRACK LENGTH
        
        data=postprocess(session_path,track_length,'yes',0);
    end
end