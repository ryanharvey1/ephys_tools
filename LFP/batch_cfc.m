% AnalyzeResults_PAE_Cylinder
clear
data=compileResults('D:\Projects\PAE_PlaceCell\ProcessedData');

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813','LEM3124'};

%% COMPILE GROUPS
data.control.id=[];
for i=1:length(control)
    data.control.id=cat(1,data.control.id,data.(control{i}).id);
end

data.pae.id=[];
for i=1:length(pae)
    data.pae.id=cat(1,data.pae.id,data.(pae{i}).id);
end

%% COMPILE GROUP IDS
group1id=data.control.id;
group2id=data.pae.id;

%GET RECORDING AREA FOR ALL SESSIONS (KUBIE) and/or TETRODES (HALO)
group1id=get_region_id(group1id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');
group2id=get_region_id(group2id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');

combined_ID=[group1id; group2id];

ca1id = unique(combined_ID(strcmp(combined_ID(:,4),'ca1'),1));
ca3id = unique(combined_ID(strcmp(combined_ID(:,4),'ca3'),1));


sess_to_run=[ca1id;ca3id];

for i=1:length(sess_to_run)
    
    disp(['running',sess_to_run{i}])
    if ~exist(['D:\Projects\PAE_PlaceCell\CFCResults\',extractBefore(sess_to_run{i},'.'),'.mat'],'file')
        try
            cfc=run_CFC(sess_to_run{i});
        catch
            disp(['running',sess_to_run{i}, ' FAILED'])
            
        end
    end
end

