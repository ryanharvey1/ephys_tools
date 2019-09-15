% populate_mazetype_helper
cd F:\ClarkP30_Recordings\AnimalMetaData
animals=dir('*.mat');
animals={animals.name};
Condition=[];

for a=1:length(animals)
    load(animals{a})
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    s=1;
    for i=size(MazeTypes,1)+1:length(sessions)+size(MazeTypes,1)
        Condition{i,1}=animals{a};
        Condition{i,2}=[sessions{s}];
        Condition{i,3}=[AnimalMetadata.RecordingLogs.(sessions{s}).MazeTypes];
        s=s+1;
    end
end
Condition
%% manually fix all maze type errors and fill everything in before continuing to next section

%% Save everything 
s=1;
for a=1:length(animals)
    load(animals{a})
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    for i=1:length(sessions)
        AnimalMetadata.RecordingLogs.(sessions{i}).ConditionTypes=Condition{s,3};
        s=s+1;
    end
    metapath=fullfile(pwd,[AnimalMetadata.AnimalName,'_metadata.mat']);
    save(metapath,'AnimalMetadata')
end



