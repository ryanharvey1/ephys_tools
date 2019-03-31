% AnalyzeResults_PAE_Cylinder_landmark_control

clear
data=compileResults('D:\Projects\PAE_PlaceCell\ProcessedData');

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813'};

%% COMPILE GROUPS
data.control.measures=[];
data.control.id=[];
for i=1:length(control)
    data.control.measures=cat(1,data.control.measures,data.(control{i}).measures);
    data.control.id=cat(1,data.control.id,data.(control{i}).id);
end

data.pae.measures=[];
data.pae.id=[];
for i=1:length(pae)
    data.pae.measures=cat(1,data.pae.measures,data.(pae{i}).measures);
    data.pae.id=cat(1,data.pae.id,data.(pae{i}).id);
end

%% COMPILE Cylinder data
session=3:4;
group1=data.control.measures(:,:,session);
group2=data.pae.measures(:,:,session);
group1id=data.control.id;
group2id=data.pae.id;


%% Delete rows with all nans (sessions without cylinder)
% to_delete=sum(isnan(group1),2)==size(group1,2);
% group1(to_delete,:,:)=[];
% group1id(to_delete,:,:)=[];
% 
% to_delete=sum(isnan(group2),2)==size(group2,2);
% group2(to_delete,:,:)=[];
% group2id(to_delete,:,:)=[];

%%
varnames=data.varnames;

%% SPLIT BY REGION
% load metadata files and extract region info
cd D:\Projects\PAE_PlaceCell\AnimalMetadata

rats=dir('*.mat');
rats={rats.name};
sess_region=[];
sessionid=[];
% mainpath='D:\Projects\PAE_PlaceCell\ProcessedData\';
mainpath=[];

for i=1:length(rats)
  load(rats{i})
  sess=fieldnames(AnimalMetadata.RecordingLogs);
  for s=1:length(sess)
      sess_region=[sess_region;{AnimalMetadata.AnimalName,sess{s},AnimalMetadata.RecordingLogs.(sess{s}).RecordingArea}];
      sessionid=[sessionid;{[mainpath,AnimalMetadata.AnimalName,'_',sess{s}]}];
  end
end

% create idx
ca1idx=strcmp(sess_region(:,3), 'ca1');
ca3idx=strcmp(sess_region(:,3), 'ca3');

ca1=sessionid(ca1idx);
ca3=sessionid(ca3idx);

% split groups between regions
% ca1
group1ca1 = group1(ismember(erase(group1id(:,1),'.mat'), ca1),:,:);
group2ca1 = group2(ismember(erase(group2id(:,1),'.mat'), ca1),:,:);
group1ca1id = group1id(ismember(erase(group1id(:,1),'.mat'), ca1),:,:);
group2ca1id = group2id(ismember(erase(group2id(:,1),'.mat'), ca1),:,:);

% ca3
group1ca3 = group1(ismember(erase(group1id(:,1),'.mat'), ca3),:,:);
group2ca3 = group2(ismember(erase(group2id(:,1),'.mat'), ca3),:,:);
group1ca3id = group1id(ismember(erase(group1id(:,1),'.mat'), ca3),:,:);
group2ca3id = group2id(ismember(erase(group2id(:,1),'.mat'), ca3),:,:);



[group1ca1,group1ca1id]=placefieldfilter(group1ca1,group1ca1id,varnames);
[group2ca1,group2ca1id]=placefieldfilter(group2ca1,group2ca1id,varnames);
[group1ca3,group1ca3id]=placefieldfilter(group1ca3,group1ca3id,varnames);
[group2ca3,group2ca3id]=placefieldfilter(group2ca3,group2ca3id,varnames);

var='Displacement';
AllStatsca1=CDFplots(group1ca1(:,ismember(varnames,var),2),...
    group2ca1(:,ismember(varnames,var),2),{'Sacc','PAE'},varnames{ismember(varnames,var)},2)

rad2deg(circ_mean(deg2rad(group1ca1(:,ismember(varnames,var),2))))
rad2deg(circ_mean(deg2rad(group2ca1(:,ismember(varnames,var),2))))

figure
polarhistogram(deg2rad(group1ca1(:,ismember(varnames,var),2)),30)

figure
polarhistogram(deg2rad(group2ca1(:,ismember(varnames,var),2)),30)

[pval table] = circ_wwtest(deg2rad(group1ca1(:,ismember(varnames,var),2)),deg2rad(group2ca1(:,ismember(varnames,var),2)))





%% Local functions below
function [group,groupid]=placefieldfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz, 
% 2) Minimum field width of 8 cm, 
% 3) Maximum field width of 80 cm, 
% 4) at least 10 trials with consistent behavior. 
% 5) at least 100 spikes

groupid=groupid(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=8 &...
    group(:,contains(varnames,'FieldWidth'))<=40 &...
    group(:,contains(varnames,'InformationContent'))>=.4 &...
    group(:,contains(varnames,'nSpikes'))>=100,:,:);

group=group(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=8 &...
    group(:,contains(varnames,'FieldWidth'))<=40 &...
    group(:,contains(varnames,'InformationContent'))>=.4 &...
    group(:,contains(varnames,'nSpikes'))>=100,:,:);
end

