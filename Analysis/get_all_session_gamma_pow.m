% get_all_session_gamma_pow
clear

path='D:\Projects\PAE_PlaceCell\ProcessedData\';

sessions=dir([path,'*.mat']);
sessions={sessions.name}';


path='D:\Projects\PAE_PlaceCell\AnimalMetadata\';
rats=dir([path,'*.mat']);
rats={rats.name}';
sessions_w_tet=[];
for i=1:length(rats)
    load(rats{i})
    disp(rats{i})
    
    current_sess=sessions(contains(sessions,extractBefore(rats{i},'_')));
    nch=AnimalMetadata.ExtracellEphys.Probes.nchannels;

    
    for s=1:length(current_sess)
        session=extractBetween(current_sess,'_','.mat');
        area=AnimalMetadata.RecordingLogs.(session{s}).RecordingArea;
        if ischar(area) || isempty(area)
            sessions_w_tet=[sessions_w_tet;[repmat(current_sess(s),nch,1),num2cell(1:nch)',repmat({area},nch,1)]];
        else
            nch=length(AnimalMetadata.RecordingLogs.(session{s}).RecordingArea);
            sessions_w_tet=[sessions_w_tet;[repmat(current_sess(s),nch,1),num2cell(1:nch)',area']];
        end
    end
    

end
combined_ID=sessions_w_tet;

% combined_ID(strcmp(combined_ID(:,3),'cortex'),:)=[];

mean_lg_freq=[];
mean_lg_pow=[];
mean_hg_freq=[];
mean_hg_pow=[];
all_area=[];
sessions=unique(combined_ID(:,1));


for i=1:length(sessions)
    load(sessions{i},'lfp');
    
%     load([extractBefore(sessions{i},'_'),'_metadata.mat'])
%     AnimalMetadata.RecordingLogs.(sessionID)
    
    disp(sessions{i})
    
    
    area=combined_ID(contains(combined_ID(:,1),sessions{i}),:);
    if size(area,1)~=size(lfp.signal,1)
        continue
    end
    
    all_area=[all_area;area];
    
    [lg_freq,lg_pow]=meanfreq(lfp.signal',1000,[25 55]);
    [hg_freq,hg_pow]=meanfreq(lfp.signal',1000,[65 120]);
    
    
    mean_lg_freq=[lg_freq';mean_lg_freq];
    mean_lg_pow=[lg_pow';mean_lg_pow];
    mean_hg_freq=[hg_freq';mean_hg_freq];
    mean_hg_pow=[hg_pow';mean_hg_pow];
end


control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813','LEM3124'};

figure
CDFplots(mean_lg_pow(contains(all_area(:,1),control) & strcmp(all_area(:,3),'ca1')),...
    mean_lg_pow(contains(all_area(:,1),pae) & strcmp(all_area(:,3),'ca1')),{'Sacc','PAE'},{'ca1 LG Power'},2)
figure
CDFplots(mean_lg_pow(contains(all_area(:,1),control) & strcmp(all_area(:,3),'ca3')),...
    mean_lg_pow(contains(all_area(:,1),pae) & strcmp(all_area(:,3),'ca3')),{'Sacc','PAE'},{'ca3 LG Power'},2)
figure
CDFplots(mean_hg_pow(contains(all_area(:,1),control) & strcmp(all_area(:,3),'ca1')),...
    mean_hg_pow(contains(all_area(:,1),pae) & strcmp(all_area(:,3),'ca1')),{'Sacc','PAE'},{'ca1 HG Power'},2)
figure
CDFplots(mean_hg_pow(contains(all_area(:,1),control) & strcmp(all_area(:,3),'ca3')),...
    mean_hg_pow(contains(all_area(:,1),pae) & strcmp(all_area(:,3),'ca3')),{'Sacc','PAE'},{'ca3 HG Power'},2)



