function groupid=get_region_id(groupid,metadatapath)
addpath(metadatapath)
rats=dir(fullfile(metadatapath,'*.mat'));
rats={rats.name};

for i=1:length(rats)
    load(rats{i})
    sess=fieldnames(AnimalMetadata.RecordingLogs);
    for s=1:length(sess)
        idx=contains(groupid(:,1),[extractBefore(rats{i},'_'),'_',sess{s}]);
        currentsess=groupid(idx,:);
        if iscell(AnimalMetadata.RecordingLogs.(sess{s}).RecordingArea)
            for r=1:length(AnimalMetadata.RecordingLogs.(sess{s}).RecordingArea)
                index = find(strcmp(currentsess(:,2), ['TT',num2str(r),'.mat']));
                currentsess(index,4)=repmat({AnimalMetadata.RecordingLogs.(sess{s}).RecordingArea{r}},length(index),1);
            end
        else
            currentsess(:,4)=repmat({AnimalMetadata.RecordingLogs.(sess{s}).RecordingArea},size(currentsess,1),1);
        end
        groupid(idx,4)=currentsess(:,4);
    end
end
end