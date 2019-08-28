% populate_mazetype_helper
% cd D:\Projects\HPCatn\AnimalMetadata
animals=dir('*.mat');
animals={animals.name};
MazeTypes=[];

for a=1:length(animals)
    load(animals{a})
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    s=1;
    for i=size(MazeTypes,1)+1:length(sessions)+size(MazeTypes,1)
        MazeTypes{i,1}=animals{a};
        MazeTypes{i,2}=[sessions{s}];
        MazeTypes{i,3}=[AnimalMetadata.RecordingLogs.(sessions{s}).MazeTypes];
        s=s+1;
    end
end
MazeTypes
%% manually fix all maze type errors and fill everything in before continuing to next section


%% or... if you have already postprocessed a bunch of sessions and need to 
% fill in the maze types based on our previous heuristic, run this section
cd ..
cd('ProcessedData')
for i=find(cellfun('isempty', MazeTypes(:,3)))'
    load([extractBefore(MazeTypes{i,1},'metadata'),MazeTypes{i,2}],'mazetypes');
    MazeTypes{i,3}=strjoin(mazetypes,',');
end
cd ..
cd('AnimalMetadata')

%% Save everything 
s=1;
for a=1:length(animals)
    load(animals{a})
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    for i=1:length(sessions)
        AnimalMetadata.RecordingLogs.(sessions{i}).MazeTypes=MazeTypes{s,3};
        s=s+1;
    end
    metapath=fullfile(pwd,[AnimalMetadata.AnimalName,'_metadata.mat']);
    save(metapath,'AnimalMetadata')
end



