% populate_mazetype_helper
% cd F:\ClarkP30_Recordings\AnimalMetaData
animals=dir('*.mat');
animals={animals.name};
Condition=[];

conditions=get_conditions;

for a=1:length(animals)
    load(animals{a})
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    s=1;
    for i=size(Condition,1)+1:length(sessions)+size(Condition,1)
        Condition{i,1}=animals{a};
        Condition{i,2}=[sessions{s}];
        Condition{i,3}=[AnimalMetadata.RecordingLogs.(sessions{s}).MazeTypes];
        Condition{i,4}=[AnimalMetadata.RecordingLogs.(sessions{s}).Notes];
        s=s+1;
    end
end


% Condition_table = table(Condition(:,1),Condition(:,2),Condition(:,3),'VariableNames',{'animal','session','maze'});
% for i=1:length(conditions)
%     Condition_table.(conditions{i})(1:size(Condition,1),1)=NaN(size(Condition,1),1);
% end



%% manually fix all maze type errors and fill everything in before continuing to next section

%% Save everything 
s=1;
for a=1:length(animals)
    load(animals{a})
    sessions=fieldnames(AnimalMetadata.RecordingLogs);
    for i=1:length(sessions)
        AnimalMetadata.RecordingLogs.(sessions{i}).ConditionTypes=Condition{s,4};
        s=s+1;
    end
    metapath=fullfile(pwd,[AnimalMetadata.AnimalName,'_metadata.mat']);
    save(metapath,'AnimalMetadata')
end



