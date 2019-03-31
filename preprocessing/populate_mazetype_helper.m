% populate_mazetype_helper
cd D:\Projects\HPCatn\AnimalMetadata
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

%%
s=1;
for a=1:length(animals)
    load(animals{a})
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    for i=1:length(sessions)
        AnimalMetadata.RecordingLogs.(sessions{i}).MazeTypes=MazeTypes{s,3};
        s=s+1;
    end
    metapath=fullfile('D:\Projects\HPCatn\AnimalMetadata',[AnimalMetadata.AnimalName,'_metadata.mat']);
    save(metapath,'AnimalMetadata')
end



