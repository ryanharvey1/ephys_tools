% AnalyzeResults_HPCatn_box
data=compileResults('D:\Projects\HPCatn\ProcessedData');

control={'HPCatn02','HPCatn05'};
lesion={'HPCatn03','HPCatn04'};

% load inactivation data
load('D:\Projects\HPCatn\AnimalMetadata\HPCatn06_metadata.mat')
sessions=fieldnames(AnimalMetadata.RecordingLogs);
for i=1:length(sessions)
    mazes{i}=AnimalMetadata.RecordingLogs.(sessions{i}).MazeTypes;
end
idx=contains(data.HPCatn06.id(:,1),sessions(contains(mazes,'box')));
tempdata=data.HPCatn06.measures(idx,:,:);
tempid=data.HPCatn06.id(idx,:);


%% COMPILE GROUPS
data.control.measures=[];
data.control.id=[];
for i=1:length(control)
    data.control.measures=cat(1,data.control.measures,data.(control{i}).measures);
    data.control.id=cat(1,data.control.id,data.(control{i}).id);
end

data.lesion.measures=[];
data.lesion.id=[];
for i=1:length(lesion)
    data.lesion.measures=cat(1,data.lesion.measures,data.(lesion{i}).measures);
    data.lesion.id=cat(1,data.lesion.id,data.(lesion{i}).id);
end
%% COMPILE box DATA
% also add index for running direction
group1=data.control.measures(:,:,3);
group2=data.lesion.measures(:,:,3);

group1id=data.control.id;
group2id=data.lesion.id;

%% add inactivation rat control data
group1=[group1;tempdata(:,:,1)];
group1id=[group1id;tempid];
control{3}='HPCatn06';

% add inactivation data
group2=[group2;tempdata(:,:,2)];
group2id=[group2id;tempid];
lesion{3}='HPCatn06';
%% DELETE MEASURES FOR OPEN ARENA
varnames=data.varnames;

colstodelete=contains(varnames,["DirectionalityIndex","Displacement","Cluster Grade",...
    "DisplacementCorr","Tightness","Incompleteness","StationInTime",...
    "TempMatch","BDistanceClust","BDistanceSpike","nlaps",...
    "rateoverlap","fieldoverlap","lap_perm_stability","stabilityoverlaps",...
    "meanstability","spatialcorrelation"]);
varnames(colstodelete)=[];
group1(:,colstodelete)=[];
group2(:,colstodelete)=[];

%% SPLIT BY REGION
% load metadata files and extract region info
cd D:\Projects\HPCatn\AnimalMetadata

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
ca3idx=contains(sess_region(:,3), ["ca3","dg"]);
cortexidx=strcmp(sess_region(:,3), 'cortex');

ca1=sessionid(ca1idx);
ca3=sessionid(ca3idx);
cortex=sessionid(cortexidx);

% split groups between regions
% ca1
group1ca1 = group1(ismember(erase(group1id(:,1),'.mat'), ca1),:,:);
group2ca1 = group2(ismember(erase(group2id(:,1),'.mat'), ca1),:,:);
group1ca1id = group1id(ismember(erase(group1id(:,1),'.mat'), ca1),:,:);
group2ca1id = group2id(ismember(erase(group2id(:,1),'.mat'), ca1),:,:);

[uCA,~,~] = uniqueRowsCA(group1ca1id);
disp([num2str(size(uCA,1)),' control ca1 cells'])
[uCA,~,~] = uniqueRowsCA(group2ca1id);
disp([num2str(size(uCA,1)),' pae ca1 cells'])

% ca3/dg
group1ca3 = group1(ismember(erase(group1id(:,1),'.mat'), ca3),:,:);
group2ca3 = group2(ismember(erase(group2id(:,1),'.mat'), ca3),:,:);
group1ca3id = group1id(ismember(erase(group1id(:,1),'.mat'), ca3),:,:);
group2ca3id = group2id(ismember(erase(group2id(:,1),'.mat'), ca3),:,:);

[uCA,~,~] = uniqueRowsCA(group1ca3id);
disp([num2str(size(uCA,1)),' control ca3 cells'])
[uCA,~,~] = uniqueRowsCA(group2ca3id);
disp([num2str(size(uCA,1)),' pae ca3 cells'])

% cortex
group1cortex = group1(ismember(erase(group1id(:,1),'.mat'), cortex),:,:);
group2cortex = group2(ismember(erase(group2id(:,1),'.mat'), cortex),:,:);
group1cortexid = group1id(ismember(erase(group1id(:,1),'.mat'), cortex),:,:);
group2cortexid = group2id(ismember(erase(group2id(:,1),'.mat'), cortex),:,:);

[uCA,~,~] = uniqueRowsCA(group1cortexid);
disp([num2str(size(uCA,1)),' control cortex cells'])
[uCA,~,~] = uniqueRowsCA(group2cortexid);
disp([num2str(size(uCA,1)),' pae cortex cells'])


%% PLACE CELL FILTER
[group1ca1,group1ca1id]=placefieldfilter(group1ca1,group1ca1id,varnames);
[group2ca1,group2ca1id]=placefieldfilter(group2ca1,group2ca1id,varnames);
[group1ca3,group1ca3id]=placefieldfilter(group1ca3,group1ca3id,varnames);
[group2ca3,group2ca3id]=placefieldfilter(group2ca3,group2ca3id,varnames);
[group1cortex,group1cortexid]=placefieldfilter(group1cortex,group1cortexid,varnames);
[group2cortex,group2cortexid]=placefieldfilter(group2cortex,group2cortexid,varnames);


[uCA,~,~] = uniqueRowsCA(group1ca1id);
disp([num2str(size(uCA,1)),' control ca1 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca1id);
disp([num2str(size(uCA,1)),' pae ca1 place cells'])

[uCA,~,~] = uniqueRowsCA(group1ca3id);
disp([num2str(size(uCA,1)),' control ca3 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca3id);
disp([num2str(size(uCA,1)),' pae ca3 place cells'])


disp('...')
disp('control rats with ca1')
C=unique(extractBefore(unique(group1ca1id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group1ca1id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('lesion rats with ca1')
C=unique(extractBefore(unique(group2ca1id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group2ca1id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('control rats with ca3')
C=unique(extractBefore(unique(group1ca3id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group1ca3id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('lesion rats with ca3')
C=unique(extractBefore(unique(group2ca3id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group2ca3id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end


% visualize_cells(group1ca1id,'D:\Projects\HPCatn\figures\placecellbox_group1_ca1')
% visualize_cells(group2ca1id,'D:\Projects\HPCatn\figures\placecellbox_group2_ca1')
% visualize_cells(group1ca3id,'D:\Projects\HPCatn\figures\placecellbox_group1_ca3')
% visualize_cells(group2ca3id,'D:\Projects\HPCatn\figures\placecellbox_group2_ca3')


group2ca1_shuff_pass=shuff(group2ca1id,'feature',{'ic','mvl','dic'},'session',2);
visualize_cells(group2ca1id(group2ca1_shuff_pass(:,1)==1,:),'D:\Projects\HPCatn\figures\placecellbox_group2_ca1_infoshuffpass')

figure
AllStatsca1=stat_plot(group1ca1,group2ca1,{'control','lesion'},varnames)
figure
AllStatsca3=stat_plot(group1ca3,group2ca3,{'control','lesion'},varnames)




%% plot individual rats
clear colors
for i=1:length(control)
    colors(:,i)=contains(group1ca1id(:,1),control{i})*i;
end
colors=sum(colors,2);

figure
plot(sort(group1ca1(colors==1,1,1)),linspace(0,1,length(group1ca1(colors==1,1,1))))
hold on
plot(sort(group1ca1(colors==2,1,1)),linspace(0,1,length(group1ca1(colors==2,1,1))))
plot(sort(group1ca1(colors==3,1,1)),linspace(0,1,length(group1ca1(colors==3,1,1))))

clear colors
for i=1:length(lesion)
    colors(:,i)=contains(group2ca1id(:,1),lesion{i})*i;
end
colors=sum(colors,2);

plot(sort(group2ca1(colors==1,1,1)),linspace(0,1,length(group2ca1(colors==1,1,1))))
plot(sort(group2ca1(colors==2,1,1)),linspace(0,1,length(group2ca1(colors==2,1,1))))
plot(sort(group2ca1(colors==3,1,1)),linspace(0,1,length(group2ca1(colors==3,1,1))))

xlabel('info content')
legend(control{1},control{2},control{3},lesion{1},lesion{2},lesion{3})

%% plot all measures to ppt
% for i=1:length(varnames)
%     AllStats=CDFplots(group1(:,i),group2(:,i),{'Control','ATN Lesion'},varnames(i),2);
%     toPPT(figure(1),'exportMode','matlab');
%     toPPT('setTitle',AllStats);
%     close all
% end

%% Locate Sessions for decoding
% clear data
% uniqueids=uniqueRowsCA(group1id);
% 
% uniquesessions=unique(group1id(:,1));
% for i=1:length(uniquesessions)
%     load(uniquesessions{i,1},'spikesID','linear_track')
%     
%     ensemble_rows=find(contains(uniqueids(:,1),uniquesessions{i,1}));
% 
%     pos=[linear_track.right{1,1}.dataspks(linear_track.right{1,1}.dataspks(:,6)==0,:);...
%         linear_track.left{1,1}.dataspks(linear_track.left{1,1}.dataspks(:,6)==0,:)];
%     [~,I]=sort(pos(:,1));
%     pos=pos(I,:);
%     
%     for k=1:length(ensemble_rows)
%         cells=find(contains(spikesID.TetrodeNum,uniqueids{ensemble_rows(k),2})...
%             & ismember(spikesID.CellNum,str2double(uniqueids(ensemble_rows(k),3))))';
%         spike_times{k,1}=sort([linear_track.right{1,cells}.dataspks(linear_track.right{1,cells}.dataspks(:,6)==1,1);...
%             linear_track.left{1,cells}.dataspks(linear_track.left{1,cells}.dataspks(:,6)==1,1)]);
%     end
%     pos_times=pos(:,1);
%     vel=pos(:,5);
%     head_angle=pos(:,4);
%     pos(:,[1,4:end])=[];
%     
%     
%     save(['D:\Projects\HPCatn\DecodingData\',uniquesessions{i,1}],'spike_times','pos','pos_times','vel','head_angle')
%     clear spike_times pos pos_times vel head_angle
% end


%%
function [group,groupid]=placefieldfilter(group,groupid,varnames)

groupid=groupid(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=9 &...
    group(:,contains(varnames,'FieldWidth'))<=40 &...
    group(:,contains(varnames,'InformationContent'))>=.3 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);

group=group(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=9 &...
    group(:,contains(varnames,'FieldWidth'))<=40 &...
    group(:,contains(varnames,'InformationContent'))>=.3 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);
end



