% function [ timestamps,StartofRec,EndofRec ] = EventSplit( path )
%
% USAGE
%
% [ timestamps,StartofRec,EndofRec ] = EventSplit( path )
%
% path              PATH TO FOLDER WITH EVENT FILE
%
%    =========================================================================
%
% OUTPUT
%
% timestamps        TIMESTAMPS OF EACH EVENT (AUTOSTART / AUTOSTOP & USER GEN: start / end)
% StartofRec        TIMESTAMPS OF START OF EACH SESSION
% EndofRec          TIMESTAMPS OF ENDS OF EACH SESSION
%
% NOTE: USER GENERATED INPUT MUST BE START & END OR START AND STOP.
% OTHERWISE EVENT TIMES WILL BE TAKEN FROM AUTOGENERATED START AND STOP TIMES
%
% RYAN HARVEY 2017

function [ timestamps,StartofRec,EndofRec ] = EventSplit( path )

% I included an option to create your own event timestamps in a mat file
% called 'events'. If you create a mat file with the above output variables
% (timestamps,StartofRec,EndofRec), named 'events.mat', that will be loaded
% here and the function will return. This is for when you mess up while
% labeling events when recording. 
% save([path filesep 'events.mat'],'timestamps','StartofRec','EndofRec')

% for Open Ephys with python script acquisition 
if ~isempty( dir([path,filesep,'**\*events_ts.csv']))
    file = dir([path,filesep,'**\*events_ts.csv']);
    events = readtable([file.folder,filesep,file.name]);
    events = events(events.Var1 >= 1,:); %keep first input thru end cause zero is null 
    StartofRec = events.start'*10^6;
    EndofRec = events.xEnd'*10^6;
    timestamps = [];
    save([path filesep 'events.mat'],'StartofRec','EndofRec','timestamps')
    return
end

if exist([path filesep 'events.mat'],'file')==2
    load([path filesep 'events.mat'])
    return
end


timestamps=[];

% LOAD EVENT FILES IF RELEVANT
year=strsplit(path,filesep); year=strsplit(char(year(end)),'-');
if str2double(year{1})>=2017; SplitPath=strsplit(path,'-');
    if str2double(SplitPath(2))>=2 || str2double(year{1})>=2017 %IF MONTH IS 2 OR MORE OR IT'S PAST 2017
        if ismac==1; [timestamps,strings,~]=Nlx2MatEV_v3([path filesep 'Events.nev'], [1 0 0 0 1], 1, 1);
        elseif ismac==0;
            try
                [timestamps,strings,~]=Nlx2MatEV([path filesep 'Events.nev'], [1 0 0 0 1], 1, 1);
            catch
                [timestamps,strings,~]=Nlx2MatEV_v4([path filesep 'Events.nev'], [1 0 0 0 1], 1, 1);
            end
        end
    end
end

% if sum(contains(strings,'Starting Recording'))==1 &&...
%         sum(contains(strings,'Stopping Recording'))==1 &&...
%         sum(contains(strings,'start'))==0 &&...
%         sum(contains(strings,'end'))==0
%     StartofRec=timestamps(1);
%     EndofRec=timestamps(2);
%     return
% end

% FIND START AND STOP TIME FROM EVENT FILES
StartofRec=1; EndofRec=1; % start will be re-written below if event exists
if isempty(timestamps)==0
    StartofRecAuto=(timestamps(ismember(strings,'Starting Recording')));
    EndofRecAuto=(timestamps(ismember(strings,'Stopping Recording')));
    if sum(ismember(strings,'start'))==0
        warning('START AND STOP TIMES NOT INDICATED IN EVENTS FILE');
        StartofRec=StartofRecAuto; EndofRec=EndofRecAuto;
    elseif sum(ismember(strings,'start'))>=1; StartofRec=(timestamps(ismember(strings,'start')'));
        if sum(ismember(strings,'end'))>=1; EndofRec=(timestamps(ismember(strings,'end')'));
        elseif sum(ismember(strings,'stop'))>=1; EndofRec=(timestamps(ismember(strings,'stop')'));
        end
    end
    
    % CHECK TO SEE IF START & END IS THE SAME LENGTH
    % IF SO, LOCATE STARTS AND ENDS OF SESSIONS DETERMINED BY INTER-TRIAL LATENCIES
    if length(StartofRec)~=length(EndofRec) && sum(ismember(strings,'start'))>=1
        warning('DIFFERENT NUMBER OF START AND END EVENTS!'); storei=1;
       for i=1:length(strings)
            if ismember(strings(i),'Stopping Recording')==1 && ismember(strings(i+1),'Starting Recording')==1
                Diff=diff(timestamps(1,i:i+1)); store(storei,:)=[Diff,i]; storei=storei+1;
            end
            if i==length(strings)-1; break; end
        end
        if exist('store','var') ==1
        [~,InterTrialIndex]=sort(store(:,1),'descend');
        try
        EndofSession=sort(timestamps(store(InterTrialIndex(1:2,1),2)'));
        catch
            StartofRec=StartofRecAuto;
            EndofRec=EndofRecAuto;
            return
        end
        StartofRec=[timestamps(1),EndofSession(1),EndofSession(2)];
        EndofRec=[EndofSession(1),EndofSession(2),timestamps(end)];
        return
        end
    end
    
    % IF AUTOSTART TO START > START TO END ---MAKE AUTOSTART THE ACTUAL START ****(DOES NOT WORK AS WRITTEN WITH UNPLUGS )****
    if ~exist('storei','var') && length(StartofRecAuto)==length(StartofRec)
        for eventnum=1:length(StartofRec)
            if StartofRec(eventnum)-StartofRecAuto(eventnum)>EndofRec(eventnum)-StartofRec(eventnum)
                disp('EVENT "START" NOT INITIATED SOON ENOUGH OR RAT GOT UNPLUGGED');
                YouGuysSureMessedUp_Didnt_You=true;
                break
            end
        end
    elseif ~exist('storei','var')
        % ALL ELSE HAS FAILED
        StartofRec=timestamps(ismember(strings,'start'));
        EndofRec=timestamps(ismember(strings,'end'));
    end
    if exist('YouGuysSureMessedUp_Didnt_You','var')==1; storei=1;
        for i=1:length(strings)
            if ismember(strings(i),'Stopping Recording')==1 && ismember(strings(i+1),'Starting Recording')==1
                Diff=diff(timestamps(1,i:i+1)); store(storei,:)=[Diff,i]; storei=storei+1;
            end
            if i==length(strings)-1; break; end
        end
        [~,InterTrialIndex]=sort(store(:,1),'descend');
        EndofSession=sort(timestamps(store(InterTrialIndex(1:2,1),2)'));
        StartofRec=[timestamps(1),EndofSession(1),EndofSession(2)];
        EndofRec=[EndofSession(1),EndofSession(2),timestamps(end)];
    end
end
end

