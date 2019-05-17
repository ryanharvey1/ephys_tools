% AnalyzeResults_PAE_LinearTrack
clear
data=compileResults('D:\Projects\PAE_PlaceCell\ProcessedData');
addpath('D:\Projects\PAE_PlaceCell\ProcessedData')

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
%% COMPILE LINEAR TRACK DATA
group1=[[data.control.measures(:,:,1),ones(size(data.control.measures(:,:,1),1),1)];...
    [data.control.measures(:,:,2),ones(size(data.control.measures(:,:,1),1),1)+1]];

group2=[[data.pae.measures(:,:,1),ones(size(data.pae.measures(:,:,1),1),1)];...
    [data.pae.measures(:,:,2),ones(size(data.pae.measures(:,:,1),1),1)+1]];

group1id=[data.control.id;data.control.id];
group2id=[data.pae.id;data.pae.id];
%% DELETE MEASURES FOR OPEN ARENA
varnames=data.varnames;
varnames=[varnames,'runningdir'];

colstodelete=contains(varnames,["Cluster Grade","borderScore","E",...
    "DisplacementCorr","bordermod","egomod","Tightness","Incompleteness",...
    "StationInTime","TempMatch","BDistanceClust","BDistanceSpike"]);
varnames(colstodelete)=[];
group1(:,colstodelete)=[];
group2(:,colstodelete)=[];


%% PLOT ALL XY FOR DEBUGING
% groupid=unique(group2id(:,1));
% for i=1:length(groupid)
%     data=load(groupid{i},'frames','events','rat','sessionID');
%     fig=figure('Name',[data.rat,'  ',data.sessionID]);
%     fig.Color=[1 1 1];
%     set(gcf, 'Position', get(0, 'Screensize'));
%     sess=size(data.events,2);
%     for ii=1:sess
%         subplot(1,sess,ii)
%         plot(data.frames(data.frames(:,1)>=data.events(1,ii) & data.frames(:,1)<=data.events(2,ii),2),...
%             data.frames(data.frames(:,1)>=data.events(1,ii) & data.frames(:,1)<=data.events(2,ii),3),'.k');
%         axis image
%     end
%     print(gcf,'-dpng', '-r80',...
%         ['D:\Projects\PAE_PlaceCell\CleanupXY',filesep,groupid{i},'.png'])
%     close all
% end

%% SPLIT BY REGION
% load metadata files and extract region info

group1id=get_region_id(group1id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');
group2id=get_region_id(group2id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');

group1ca1 = group1(strcmp(group1id(:,4),'ca1'),:);
group1ca1id = group1id(strcmp(group1id(:,4),'ca1'),:);
group1ca3 = group1(strcmp(group1id(:,4),'ca3'),:);
group1ca3id = group1id(strcmp(group1id(:,4),'ca3'),:);
group1cortex = group1(strcmp(group1id(:,4),'cortex'),:);
group1cortexid = group1id(strcmp(group1id(:,4),'cortex'),:);

group2ca1 = group2(strcmp(group2id(:,4),'ca1'),:);
group2ca1id = group2id(strcmp(group2id(:,4),'ca1'),:);
group2ca3 = group2(strcmp(group2id(:,4),'ca3'),:);
group2ca3id = group2id(strcmp(group2id(:,4),'ca3'),:);
group2cortex = group2(strcmp(group2id(:,4),'cortex'),:);
group2cortexid = group2id(strcmp(group2id(:,4),'cortex'),:);


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

%% cluster dentate from CA3
% x=group1ca3(:,:,1);
% x(any(isnan(x) | isinf(x),2),:)=[];
% x=zscore(x);
% [coeff,score,latent,tsquared,explained,mu]=pca(x);
% 
% figure
% scatter3(score(:,1),score(:,2),score(:,3))
% axis equal
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')
% 
% scatter(group1ca3(:,contains(varnames,'spikewidth'),1),group1ca3(:,contains(varnames,'burstIdx'),1),'k','Filled')
% 


%% PLACE CELL FILTER
[group1ca1,group1ca1id]=placefieldfilter(group1ca1,group1ca1id,varnames);
[group2ca1,group2ca1id]=placefieldfilter(group2ca1,group2ca1id,varnames);
[group1ca3,group1ca3id]=placefieldfilter(group1ca3,group1ca3id,varnames);
[group2ca3,group2ca3id]=placefieldfilter(group2ca3,group2ca3id,varnames);


[uCA,~,~] = uniqueRowsCA(group1ca1id);
disp([num2str(size(uCA,1)),' control ca1 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca1id);
disp([num2str(size(uCA,1)),' pae ca1 place cells'])

[uCA,~,~] = uniqueRowsCA(group1ca3id);
disp([num2str(size(uCA,1)),' control ca3 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca3id);
disp([num2str(size(uCA,1)),' pae ca3 place cells'])

uniqueRowsCA(extractBefore(group2ca3id(:,1),'_'))


%  [pontential]=COM(group1ca1id)
% visualizecells(uniqueRowsCA(group1ca1id),'sacca1')
% visualizecells(uniqueRowsCA(group2ca1id),'paeca1')
% visualizecells(uniqueRowsCA(group1ca3id),'sacca3')
% visualizecells(uniqueRowsCA(group2ca3id),'paeca3')

fig=figure;fig.Color=[1 1 1];
AllStatsca3=stat_plot(group1ca3,group2ca3,{'Sacc','PAE'},varnames,'plots',1,'plottype','beeswarm')

% AllStatsca1=CDFplots(group1ca1,group2ca1,{'Sacc','PAE'},varnames,1)
    AllStatsca3=CDFplots(group1ca3,group2ca3,{'Sacc','PAE'},varnames,1);
    
    ph1=group1ca3(group1ca3(:,contains(varnames,'Phpval'))<.05,contains(varnames,'PhcircLinCorr'));
    ph2=group2ca3(group2ca3(:,contains(varnames,'Phpval'))<.05,contains(varnames,'PhcircLinCorr'));
    figure
    PHca3=CDFplots(ph1,ph2,{'Sacc','PAE'},varnames{contains(varnames,'PhcircLinCorr')},1)
fig=figure;fig.Color=[1 1 1]
[h1,h2]=CoolHistogram(ph1,ph2,50,varnames{contains(varnames,'PhcircLinCorr')})


for i=1:length(varnames)
    fig=figure('Name',['ca3 ',varnames{i}],'NumberTitle','off');
    AllStatsca3=CDFplots(group1ca3(:,i),group2ca3(:,i),{'Sacc','PAE'},varnames{i},2);
    toPPT(fig,'exportMode','matlab');
    toPPT('setTitle',AllStatsca3);
    close all
end

visualize_cells(uniqueRowsCA(group1ca1id),...
    'D:\Projects\PAE_PlaceCell\Figures\PlaceCells\control_ca1');
visualize_cells(uniqueRowsCA(group2ca1id),...
    'D:\Projects\PAE_PlaceCell\Figures\PlaceCells\pae_ca1');
visualize_cells(uniqueRowsCA(group1ca3id),...
    'D:\Projects\PAE_PlaceCell\Figures\PlaceCells\control_ca3');
visualize_cells(uniqueRowsCA(group2ca3id),...
    'D:\Projects\PAE_PlaceCell\Figures\PlaceCells\pae_ca3');

% visualizecells(group1ca1id,'control_ca1')
% visualizecells(group2ca1id,'pae_ca1')
% visualizecells(group1ca3id,'control_ca3')
% visualizecells(group2ca3id,'pae_ca3')


popvector(group1ca1id,group1ca1(:,end),'controlca1');
popvector(group1ca3id,group1ca3(:,end),'controlca3');

popvector(group2ca1id,group2ca1(:,end),'paeca1');
popvector(group2ca3id,group2ca3(:,end),'paeca3');


% AllStatsca1=CDFplots(group1ca1(:,contains(varnames,'thetaindex')),group2ca1(:,contains(varnames,'thetaindex')),{'Sacc','PAE'},'thetaindex',1)



% close all
% [pontential]=findexamples(group2ca1,group2ca1id,varnames,'thetaindex');
% close all
% [pontential]=findexamples(group1ca1,group1ca1id,varnames,'thetaindex');


% print(gcf,'-dmeta',['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\PAE Project\Presentations\SfN2018\CellExamples',...
%     filesep,'LS19_S20170508130713.emf'])

% data=load('LS19_S20170508130713.mat')
% postprocessFigures(data,{'TT4.mat',8});
% 
% groupid= group1ca1id;
% for i=1:length(groupid)
%     concatid{i,1}=strcat(groupid{i,:});
% end
% [~,~,ic]=unique(concatid,'stable');
% idx=find(abs(diff(ic))>1);
% idx=idx(2);
% 
% for i=1:length(group1ca1id)
% if i>idx
%     d=2;
% else
%     d=1;
% end
% load(group1ca1id{i,1},'thetaautocorr','spikesID')
% cells=find(contains(spikesID.TetrodeNum,group1ca1id{i,2}) & ismember(spikesID.CellNum,str2double(group1ca1id{i,3})))';
% autocorr(i,:)=thetaautocorr{cells,d}';
% end
% group1ca1(group1ca1(:,contains(varnames,'thetaindex'))>1.5,contains(varnames,'thetaindex'))
% 
% [~,I]=sort(group1ca1(:,contains(varnames,'thetaindex')));
% figure;
% imagesc(autocorr(I,:))



% p=[];
% for i=1:1000
%     try
%     AllStatsca3=CDFplots(group1ca3(randi([1,size(group1ca3,1)],13,1),:),group2ca3,{'Sacc','PAE'},varnames,3);
%     p=[ p,str2double(extractBetween(AllStatsca3,'p=',','))];
%     catch
%     end
% end
% (sum(p<.05,2)/size(p,2))>.95

AllStatsca1=CDFplots(group1ca1(:,1),group2ca1(:,1),{'Sacc','PAE'},varnames{1},2)
AllStatsca1=CDFplots(group1ca1(:,2),group2ca1(:,2),{'Sacc','PAE'},varnames{2},2)
AllStatsca1=CDFplots(group1ca1(:,4),group2ca1(:,4),{'Sacc','PAE'},varnames{4},2)
AllStatsca1=CDFplots(group1ca1(:,5),group2ca1(:,5),{'Sacc','PAE'},varnames{5},2)

AllStatsca1=CDFplots(group1ca1(:,19),group2ca1(:,19),{'Sacc','PAE'},varnames{19},2)

AllStatsca1=CDFplots(group1ca1(group1ca1(:,20)<.05,19),group2ca1(group2ca1(:,20)<.05,19),{'Sacc','PAE'},varnames{19},2)

fig=figure;fig.Color=[1 1 1]
[h1,h2]=CoolHistogram(group1ca1(group1ca1(:,20)<.05,19),group2ca1(group2ca1(:,20)<.05,19),50,varnames{19})




%%
[group1ca1_stab]=lapstability(group1ca1id);
[group2ca1_stab]=lapstability(group2ca1id);
[group1ca3_stab]=lapstability(group1ca3id);
[group2ca3_stab]=lapstability(group2ca3id);


RL_anova(group1ca1_stab.cor,group2ca1_stab.cor,'ca1 Spatial Corr')
RL_anova(group1ca3_stab.cor,group2ca3_stab.cor,'ca3 Spatial Corr')

RL_anova(group1ca1_stab.rate,group2ca1_stab.rate,'ca1 rate')
RL_anova(group1ca3_stab.rate,group2ca3_stab.rate,'ca3 rate')

RL_anova(group1ca1_stab.spatialstd,group2ca1_stab.spatialstd,'ca1 std')
RL_anova(group1ca3_stab.spatialstd,group2ca3_stab.spatialstd,'ca3 std')

RL_anova(group1ca1_stab.com,group2ca1_stab.com,'ca1 com')
RL_anova(group1ca3_stab.com,group2ca3_stab.com,'ca3 com')


CDFplots(nanstd(group1ca1_stab.com,0,2)/sqrt(size(group1ca1_stab.com,1)),...
    nanstd(group2ca1_stab.com,0,2)/sqrt(size(group2ca1_stab.com,1)),{'Sacc','PAE'},{'std_com'},2)

CDFplots(nanstd(group1ca3_stab.com,0,2)/sqrt(size(group1ca3_stab.com,1)),...
    nanstd(group2ca3_stab.com,0,2)/sqrt(size(group2ca3_stab.com,1)),{'Sacc','PAE'},{'std_com'},2)
%%
% function [pontential]=COM(groupid)
% cd D:\Projects\PAE_PlaceCell\ProcessedData
% for i=1:length(groupid)
% load(groupid{i,1},'ratemap','spikesID');
%     cells=find(contains(spikesID.TetrodeNum,groupid{i,2}) &...
%         ismember(spikesID.CellNum,str2double(groupid{i,3})))';
%     ratemap(cells,:,1)
%     ratemap(cells,:,2)
% [fields1]=getPlaceFields(ratemap{cells,:,1},'minPeakRate',0,'minFieldWidth',3,'maxFieldWidth',length(ratemap{cells,:,1}));
% [fields2]=getPlaceFields(ratemap{cells,:,2},'minPeakRate',0,'minFieldWidth',3,'maxFieldWidth',length(ratemap{cells,:,2}));
% 
%                 for fie=1:length(fields1{1, 1})
%                     PR(fie)=fields1{1, 1}{fie}.peakFR;
%                 end
% end
% end
function [pontential]=findexamples(groupdata,groupid,varnames,var)
cd D:\Projects\PAE_PlaceCell\ProcessedData
x=groupdata(:,contains(varnames,var));
meanx=mean(x);
stdx=std(x);
disp([num2str(meanx),'+-',num2str(stdx)])
% pontential=groupid(x<meanx+stdx*.5 & x>meanx-stdx*.5,:);

pontential=groupid(x>meanx+stdx*.5,:);

pontential=pontential(randi([1,size(pontential,1)],10,1),:);
for i=1:length(pontential)
    data=load(pontential{i,1});
    cells=find(contains(data.spikesID.TetrodeNum,pontential{i,2}) & ismember(data.spikesID.CellNum,str2double(pontential{i,3})))';
    disp([data.sessionID,': ',var,' ',num2str(squeeze(data.measures(cells,contains(data.varnames,var),:))')])
    postprocessFigures(data,{pontential{i,2},str2double(pontential(i,3))});
end
end

function [results]=lapstability(groupid)
cd D:\Projects\PAE_PlaceCell\ProcessedData
for i=1:length(groupid)
    concatid{i,1}=strcat(groupid{i,:});
end
[~,~,ic]=unique(concatid,'stable');
idx=find(abs(diff(ic))>1);
idx=idx(2);

for i=1:length(groupid)
    if i<idx
        direction='right';
    else
        direction='left';
    end
    load(groupid{i,1},'linear_track','spikesID');
    cells=find(contains(spikesID.TetrodeNum,groupid{i,2}) & ismember(spikesID.CellNum,str2double(groupid{i,3})))';
    laps=vertcat(linear_track.(direction){cells}.maps{:});
    if i<idx
       laps=fliplr(laps); 
    end
    laps=imresize(laps,[NaN,40]);
    laps=rescale(laps,0,1);
    overallmap=mean(laps,1);
    if size(laps,1)>=10
        [~,I]=max(laps(1:10,:),[],2);
        com(i,:)=diff(I);
        spatialstd(i,:)=nanstd(laps(1:10,:),0,1);
        for ilap=1:10
            cor(i,ilap)=corr2(laps(ilap,:),overallmap);
            rate(i,ilap)=abs(max(laps(ilap,:))-max(overallmap));
        end
    end
end
results.cor=cor;
results.rate=rate;
results.spatialstd=spatialstd;
results.com=com;

rowstodelete=sum(results.cor==0,2)==size(results.cor,2);
results.cor(rowstodelete,:)=[];
results.rate(rowstodelete,:)=[];
results.spatialstd(rowstodelete,:)=[];
results.com(rowstodelete,:)=[];
end






function [group,groupid]=placefieldfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 1 Hz, 
% 2) Minimum field width of 8 cm, 
% 3) Maximum field width of 80 cm, 
% 4) at least 10 trials with consistent behavior. 
% 5) at least 100 spikes

groupid=groupid(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=8 &...
    group(:,contains(varnames,'FieldWidth'))<=80 &...
    group(:,contains(varnames,'nlaps'))>=10 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);

group=group(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=8 &...
    group(:,contains(varnames,'FieldWidth'))<=80 &...
    group(:,contains(varnames,'nlaps'))>=10 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);
end

function visualizecells(groupid,group)
cd('D:\Projects\PAE_PlaceCell\ProcessedData')
for i=1:length(groupid)
    data=load(groupid{i,1});
    try
        postprocessFigures.main(data,{groupid{i,2},str2double(groupid(i,3))});
        pause(.001)
    catch
    end
%     set(gcf, 'Position', get(0, 'Screensize'));

    set(gcf,'WindowState','maximized')

    print(gcf,'-dpng', '-r70',...
        ['C:\Users\ryanh\Dropbox\school work\UNM\Lab\Projects\PAE Project\Presentations\SfN2019\CellExamples\',group,...
        filesep,groupid{i,1},groupid{i,2},groupid{i,3},'.png'])
    close all
end
end

function popvector(groupid,runningdir,group)
for i=1:length(groupid)
    data=load(groupid{i,1},'ratemap','spikesID');
    cells=find(contains(data.spikesID.TetrodeNum,groupid{i,2}) &...
        ismember(data.spikesID.CellNum,str2double(groupid{i,3})))';

    placecellmaps(i,:)=imresize(rescale(data.ratemap{cells,runningdir(i)},0,1),[1,40]);
    
    if runningdir(i)==1
        oppositerunningmaps(i,:)=imresize(rescale(data.ratemap{cells,2},0,1),[1,40]);
    else
        oppositerunningmaps(i,:)=imresize(rescale(data.ratemap{cells,1},0,1),[1,40]);
    end


end

[~,I]=max(placecellmaps,[],2);
[~,I2]=sort(I);

placecellmaps=placecellmaps(I2,:);
oppositerunningmaps=oppositerunningmaps(I2,:);

fig=figure;fig.Color=[1 1 1];
subplot(1,2,1)
imagesc(placecellmaps)
ylabel('cells')
xlabel('position bins')
set(gca,'FontSize',12)

subplot(1,2,2)
imagesc(oppositerunningmaps)
ylabel('cells')
xlabel('position bins')
set(gca,'FontSize',12)

set(gcf, 'Position', get(0, 'Screensize'));
print(gcf,'-dpng', '-r400',...
    ['D:\Projects\PAE_PlaceCell\Figures\popvectors\',group,filesep,'popvec.png'])

savefig(['D:\Projects\PAE_PlaceCell\Figures\popvectors\',group,filesep,'popvec.fig'])
close all
end