
% batch_cfc
cd D:\Projects\PAE_PlaceCell\AnimalMetadata

rats=dir('*.mat');
rats={rats.name};
sess_region=[];
sessionid=[];
mainpath='D:\Projects\PAE_PlaceCell\ProcessedData\';
for i=1:length(rats)
  load(rats{i})
  sess=fieldnames(AnimalMetadata.RecordingLogs);
  for s=1:length(sess)
      sess_region=[sess_region;{AnimalMetadata.AnimalName,sess{s},AnimalMetadata.RecordingLogs.(sess{s}).RecordingArea}];
      sessionid=[sessionid;{[mainpath,AnimalMetadata.AnimalName,'_',sess{s}]}];
  end
end



% sum(strcmp(sess_region(:,3), 'cortex'))
ca1idx=strcmp(sess_region(:,3), 'ca1');
ca3idx=strcmp(sess_region(:,3), 'ca3');


sess_to_run=sessionid(ca1idx | ca3idx);
for i=1:length(sess_to_run)
    disp(['running',sess_to_run{i}])
    session=strsplit(sess_to_run{i},filesep);
    if ~exist(['D:\Projects\PAE_PlaceCell\CFCResults\',session{end},'.mat'],'file')
        
        %
        
        cfc=run_CFC(sess_to_run{i});
        
    end
end

