%% Assign_region_to_metadata

% clear
% load metadata of choice
% load('D:\Projects\PAE_PlaceCell\AnimalMetadata\LEM3120_metadata.mat')

%% First, fill in recording logs based on processed data files
% cd to your processed data folder

if contains(pwd,'ProcessedData')~=1
    error('cd to ProcessedData')
end

files=dir;

files={files.name};

files=files(contains(files,AnimalMetadata.AnimalName));

files=extractBetween(files,'_','.mat');

for i=1:length(files)
    if isfield(AnimalMetadata.RecordingLogs,(files{i}))
        continue
    end
    AnimalMetadata.RecordingLogs.(files{i}).MazeTypes='';
    AnimalMetadata.RecordingLogs.(files{i}).ConditionTypes='';
    AnimalMetadata.RecordingLogs.(files{i}).DorsalVentral=[];
    AnimalMetadata.RecordingLogs.(files{i}).RecordingArea='';
    AnimalMetadata.RecordingLogs.(files{i}).Notes='';
end

%% Secondly, run this session and open up areas var (it's helpful to undock areas and half screen it)
%  fill in missing regions 
%  Importantly, You can load and plot individual sessions with the
%  following section to check waveform and other characteristics to help in
%  assigning region ids
sessions=fieldnames(AnimalMetadata.RecordingLogs);
for i=1:length(sessions)
    areas{i,1}=[sessions{i}];
    areas{i,2}=[AnimalMetadata.RecordingLogs.(sessions{i}).RecordingArea];
end

%% load and plot
% enter row number of session you want to look at below 
row_num=22;

close all
disp(['Loading...',[AnimalMetadata.AnimalName,'_',areas{row_num,1}]])
data=load([AnimalMetadata.AnimalName,'_',areas{row_num,1}]);
disp('Plotting...')
postprocessFigures.main(data);

% % plot and save all
% for i=1:length(areas)
%     disp(['Loading...',[AnimalMetadata.AnimalName,'_',areas{i,1}]])
%     data=load([AnimalMetadata.AnimalName,'_',areas{i,1}]);
%     disp('Plotting...')
%     postprocessFigures.main(data);
%     mkdir(['D:\Projects\PAE_PlaceCell\Figures\LEM3124\',data.sessionID])
% 
%     FigList=findobj(allchild(0), 'flat', 'Type', 'figure');
%     for iFig=1:length(FigList)
%         FigHandle=FigList(iFig);
%         FigName=get(FigHandle, 'Name');
%         
%         set(FigHandle, 'Position', get(0, 'Screensize'));
%         erase(FigName,{' ',':'})
%         saveas(FigHandle,['D:\Projects\PAE_PlaceCell\Figures\LEM3124\',data.sessionID,'\',erase(FigName,{' ',':'}),'.jpeg'])
%     end
%     close all
% end



%% Third, run this session when you are done filling out the mission regions
sessions=fieldnames(AnimalMetadata.RecordingLogs);
for i=1:length(sessions)
    AnimalMetadata.RecordingLogs.(sessions{i}).RecordingArea=char(areas{i,2:end});
end

%% Fourth, save the updated metadata file 
cd ..
% metapath=strsplit(dataset,filesep);
metapath = pwd;
metapath=fullfile(metapath,'AnimalMetadata',[AnimalMetadata.AnimalName,'_metadata']);
save(metapath,'AnimalMetadata')





