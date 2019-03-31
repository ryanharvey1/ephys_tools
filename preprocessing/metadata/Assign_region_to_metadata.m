%% Assign_region_to_metadata
clear
% load metadata of choice
load('D:\Projects\HPCatn\AnimalMetadata\HPCatn06_metadata.mat')

%% First, fill in recording logs based on processed data files

dataset='D:\Projects\HPCatn\ProcessedData';
cd(dataset)

files=dir;

files={files.name};

files=files(contains(files,AnimalMetadata.AnimalName));

files=extractBetween(files,'_','.mat');

for i=1:length(files)
    if isfield(AnimalMetadata.RecordingLogs,(files{i}))
        continue
    end
    AnimalMetadata.RecordingLogs.(files{i}).MazeTypes='';
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
postprocessFigures_v3(data);


%% Third, run this session when you are done filling out the mission regions
sessions=fieldnames(AnimalMetadata.RecordingLogs);
for i=1:length(sessions)
    AnimalMetadata.RecordingLogs.(sessions{i}).RecordingArea=areas{i,2};
end

%% Fourth, save the updated metadata file 
metapath=strsplit(dataset,filesep);
metapath=fullfile(metapath{1:end-1},'AnimalMetadata',[AnimalMetadata.AnimalName,'_metadata']);
save(metapath,'AnimalMetadata')





