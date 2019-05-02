% AnalyzeResults_PAE_Cylinder_landmark_control

clear
data=compileResults('D:\Projects\PAE_PlaceCell\ProcessedData');

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813','LEM3124'};

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
group1id=get_region_id(group1id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');
group2id=get_region_id(group2id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');

group1ca1 = group1(strcmp(group1id(:,4),'ca1'),:,:);
group1ca1id = group1id(strcmp(group1id(:,4),'ca1'),:,:);
group1ca3 = group1(strcmp(group1id(:,4),'ca3'),:,:);
group1ca3id = group1id(strcmp(group1id(:,4),'ca3'),:,:);
group1cortex = group1(strcmp(group1id(:,4),'cortex'),:,:);
group1cortexid = group1id(strcmp(group1id(:,4),'cortex'),:,:);

group2ca1 = group2(strcmp(group2id(:,4),'ca1'),:,:);
group2ca1id = group2id(strcmp(group2id(:,4),'ca1'),:,:);
group2ca3 = group2(strcmp(group2id(:,4),'ca3'),:,:);
group2ca3id = group2id(strcmp(group2id(:,4),'ca3'),:,:);
group2cortex = group2(strcmp(group2id(:,4),'cortex'),:,:);
group2cortexid = group2id(strcmp(group2id(:,4),'cortex'),:,:);

[uCA,~,~] = uniqueRowsCA(group1ca1id);
disp([num2str(size(uCA,1)),' control ca1 cells'])
[uCA,~,~] = uniqueRowsCA(group2ca1id);
disp([num2str(size(uCA,1)),' pae ca1 cells'])

[uCA,~,~] = uniqueRowsCA(group1ca3id);
disp([num2str(size(uCA,1)),' control ca3 cells'])
[uCA,~,~] = uniqueRowsCA(group2ca3id);
disp([num2str(size(uCA,1)),' pae ca3 cells'])

[uCA,~,~] = uniqueRowsCA(group1cortexid);
disp([num2str(size(uCA,1)),' control cortex cells'])
[uCA,~,~] = uniqueRowsCA(group2cortexid);
disp([num2str(size(uCA,1)),' pae cortex cells'])

%% place cell filter
[group1ca1,group1ca1id]=placefieldfilter(group1ca1,group1ca1id,varnames);
[group2ca1,group2ca1id]=placefieldfilter(group2ca1,group2ca1id,varnames);
[group1ca3,group1ca3id]=placefieldfilter(group1ca3,group1ca3id,varnames);
[group2ca3,group2ca3id]=placefieldfilter(group2ca3,group2ca3id,varnames);

var='Displacement';


% AllStatsca1=CDFplots(group1ca1(:,ismember(varnames,var),2),...
%     group2ca1(:,ismember(varnames,var),2),{'Sacc','PAE'},varnames{ismember(varnames,var)},2)

rad2deg(circ_mean(deg2rad(group1ca1(:,ismember(varnames,var),2))))
rad2deg(circ_mean(deg2rad(group2ca1(:,ismember(varnames,var),2))))

figure
subplot(2,2,1)
polarhistogram(deg2rad(group1ca1(group1ca1(:,ismember(varnames,'DisplacementCorr'),2)>.3,...
    ismember(varnames,var),2)),30)
title('control ca1 landmark control')

subplot(2,2,2)
polarhistogram(deg2rad(group2ca1(group2ca1(:,ismember(varnames,'DisplacementCorr'),2)>.3,...
    ismember(varnames,var),2)),30)
title('pae ca1 landmark control')

subplot(2,2,3)
polarhistogram(deg2rad(group1ca3(group1ca3(:,ismember(varnames,'DisplacementCorr'),2)>.3,...
    ismember(varnames,var),2)),30)
title('control ca3 landmark control')

subplot(2,2,4)
polarhistogram(deg2rad(group2ca3(group2ca3(:,ismember(varnames,'DisplacementCorr'),2)>.3,...
    ismember(varnames,var),2)),30)
title('pae ca3 landmark control')




[pval table] = circ_wwtest(deg2rad(group1ca1(:,ismember(varnames,var),2)),...
    deg2rad(group2ca1(:,ismember(varnames,var),2)))

% group2ca1id(28,:)
% i=8
% 
%  data=load(group2ca1id{i,1});
%     postprocessFigures.main(data,{group2ca1id{i,2},str2double(group2ca1id(i,3))});



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
    group(:,contains(varnames,'InformationContent'))>=.3 &...
    group(:,contains(varnames,'nSpikes'))>=100,:,:);

group=group(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=8 &...
    group(:,contains(varnames,'FieldWidth'))<=40 &...
    group(:,contains(varnames,'InformationContent'))>=.3 &...
    group(:,contains(varnames,'nSpikes'))>=100,:,:);
end

