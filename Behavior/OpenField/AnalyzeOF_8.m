
%%
%%%%%%%%%%% Plot Paths with home base outline and cue location

load('params_V16.mat')
A = table2array(params);

for i=1:size(A,1)
   for j=1:size(A,2) 
       classes{i,j}=class(A{i,j});
       sizes(i,j)=size(A{i,j},1)==1 & size(A{i,j},2)==1;
   end
end

idx=contains(classes,'double') & sizes;

idx=sum(idx,1)==size(idx,1);

datamat=cell2mat(A(:,idx));

vars=fieldnames(params);

variables=vars(idx);

id=split(params.subID,'_');

day=id(:,2);

[~,~,subject]=unique(id(:,1))

id=id(:,1);

%% Get Home Base critera 
numHb=cellfun('size',params.hbOcc,2); 

for i=1:size(params.HBcenter,1)
    [~,occIdx(i,1)]=max(params.hbOcc{i});
    primary_hbOcc(i,1)=params.hbOcc{i}(1,occIdx(i,1));
    primary_hbEntry(i,1)=params.entries{i}{1,occIdx(i,1)};
    if ~isempty(params.cueCoords{i})
    primary_hbDist2Cue(i,1)=params.HBdist2Cue{i}{1,occIdx(i,1)};
    else
    primary_hbDist2Cue(i,1)=NaN;
    end
    primary_hbStops(i,1)=params.HBstops{i}{1,occIdx(i,1)};
    primary_hbvel(i,1)=params.hbVel{i}(1,occIdx(i,1));
    primary_area(i,1)=params.fieldarea{i}(1,occIdx(i,1));
end

vars={'Subject','day','group','NumHBs','Occupancy','Entries','Distance2Cue','NumStops','AvgVelocity','area'};
HBmeas=[numHb primary_hbOcc primary_hbEntry primary_hbDist2Cue primary_hbStops primary_hbvel];
day=repmat({'D1';'D2'},size(params.HBcenter,1)/2,1);
group=[repmat({'Tg'},size(params.HBcenter,1)/2,1);repmat({'Wt'},size(params.HBcenter,1)/2,1)];
id=split(params.subID,'_');
[~,~,subject]=unique(id(:,1));

hbMeas2Python=table(subject,day,group,numHb,primary_hbOcc,primary_hbEntry, primary_hbDist2Cue,...
    primary_hbStops, primary_hbvel,primary_area,'VariableNames',vars);

writetable(hbMeas2Python,'D:\Users\BClarkLab\Google Drive (lberkowitz@unm.edu)\Manuscripts\In Progress\TgF344-AD_OF\Data\hbData.csv');

%% Proximity of primary hb day 1 versus primary hb day 2
row=1;
for i=1:2:size(params,1)
    primaryHBdist(row,1)=sqrt((params.HBcenter{i,1}{1,occIdx(i,1)}(1,1)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,1))^2+...
        (params.HBcenter{i,1}{1,occIdx(i,1)}(1,2)-params.HBcenter{i+1,1}{1,occIdx(i+1,1)}(1,2))^2);
    row=row+1;
end

figure; 
plotspread_wrapper(primaryHBdist(13:end,1),primaryHBdist(1:12,1),{'WT','Tg'})
title('Distance between Day 1 & Day 2 Primary Home Base Centroids')

