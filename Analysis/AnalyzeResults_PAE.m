% AnalyzeResults_PAE_LinearTrack
clear
data=compileResults('F:\Projects\PAE_PlaceCell\ProcessedData');
addpath('F:\Projects\PAE_PlaceCell\ProcessedData')

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120','LEM3216'};
pae={'RH11','RH16','LS17','LS19','LE2813','LEM3124','LEM3206','LEM3246'};

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

% colstodelete=contains(varnames,["Cluster Grade","borderScore","E",...
%     "DisplacementCorr","bordermod","egomod","Tightness","Incompleteness",...
%     "StationInTime","TempMatch","BDistanceClust","BDistanceSpike"]);
% varnames(colstodelete)=[];
% group1(:,colstodelete)=[];
% group2(:,colstodelete)=[];


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

group1id=get_region_id(group1id,'F:\Projects\PAE_PlaceCell\AnimalMetadata');
group2id=get_region_id(group2id,'F:\Projects\PAE_PlaceCell\AnimalMetadata');

group1ca1 = group1(strcmp(group1id(:,4),'ca1'),:);
group1ca1id = group1id(strcmp(group1id(:,4),'ca1'),:);
group1ca3 = group1(strcmp(group1id(:,4),'ca3'),:);
group1ca3id = group1id(strcmp(group1id(:,4),'ca3'),:);
group1dg = group1(strcmp(group1id(:,4),'dg'),:);
group1dgid = group1id(strcmp(group1id(:,4),'dg'),:);
group1cortex = group1(strcmp(group1id(:,4),'cortex'),:);
group1cortexid = group1id(strcmp(group1id(:,4),'cortex'),:);

group2ca1 = group2(strcmp(group2id(:,4),'ca1'),:);
group2ca1id = group2id(strcmp(group2id(:,4),'ca1'),:);
group2ca3 = group2(strcmp(group2id(:,4),'ca3'),:);
group2ca3id = group2id(strcmp(group2id(:,4),'ca3'),:);
group2dg = group2(strcmp(group2id(:,4),'dg'),:);
group2dgid = group2id(strcmp(group2id(:,4),'dg'),:);
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

[uCA,~,~] = uniqueRowsCA(group1dgid);
disp([num2str(size(uCA,1)),' control dg cells'])
[uCA,~,~] = uniqueRowsCA(group2dgid);
disp([num2str(size(uCA,1)),' pae dg cells'])

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
% [group1ca1,group1ca1id]=preshuffle_placefieldfilter(group1ca1,group1ca1id,varnames);
% [group2ca1,group2ca1id]=preshuffle_placefieldfilter(group2ca1,group2ca1id,varnames);
% [group1ca3,group1ca3id]=preshuffle_placefieldfilter(group1ca3,group1ca3id,varnames);
% [group2ca3,group2ca3id]=preshuffle_placefieldfilter(group2ca3,group2ca3id,varnames);
% [group1dg,group1dgid]=preshuffle_placefieldfilter(group1dg,group1dgid,varnames);
% [group2dg,group2dgid]=preshuffle_placefieldfilter(group2dg,group2dgid,varnames);
% 
% 
% group1ca1_shuff_pass=shuff(group1ca1id,'runningdir',group1ca1(:,contains(varnames,'runningdir')),'feature',{'ic'});
% group2ca1_shuff_pass=shuff(group2ca1id,'runningdir',group2ca1(:,contains(varnames,'runningdir')),'feature',{'ic'});
% group1ca3_shuff_pass=shuff(group1ca3id,'runningdir',group1ca3(:,contains(varnames,'runningdir')),'feature',{'ic'});
% group2ca3_shuff_pass=shuff(group2ca3id,'runningdir',group2ca3(:,contains(varnames,'runningdir')),'feature',{'ic'});
% group1dg_shuff_pass=shuff(group1dgid,'runningdir',group1dg(:,contains(varnames,'runningdir')),'feature',{'ic'});
% group2dg_shuff_pass=shuff(group2dgid,'runningdir',group2dg(:,contains(varnames,'runningdir')),'feature',{'ic'});
% 
% 
% fig=figure;fig.Color=[1 1 1];
% AllStatsca1=stat_plot(group1ca1(logical(group1ca1_shuff_pass),:),...
%     group2ca1(logical(group2ca1_shuff_pass),:),{'Sacc','PAE'},varnames)
% AllStatsca1=stat_plot(group1ca3(logical(group1ca3_shuff_pass),:),...
%     group2ca3(logical(group2ca3_shuff_pass),:),{'Sacc','PAE'},varnames)
% 
% visualize_cells(uniqueRowsCA(group1ca1id(logical(group1ca1_shuff_pass),:)),...
%     'D:\Projects\PAE_PlaceCell\Figures\shufflepass\group1ca1_shuff_pass','dpi','-r80')
% visualize_cells(uniqueRowsCA(group2ca1id(logical(group2ca1_shuff_pass),:)),...
%     'D:\Projects\PAE_PlaceCell\Figures\shufflepass\group2ca1_shuff_pass','dpi','-r80')
% visualize_cells(uniqueRowsCA(group1ca3id(logical(group1ca3_shuff_pass),:)),...
%     'D:\Projects\PAE_PlaceCell\Figures\shufflepass\group1ca3_shuff_pass','dpi','-r80')
% visualize_cells(uniqueRowsCA(group2ca3id(logical(group2ca3_shuff_pass),:)),...
%     'D:\Projects\PAE_PlaceCell\Figures\shufflepass\group2ca3_shuff_pass','dpi','-r80')
% visualize_cells(uniqueRowsCA(group1dgid(logical(group1dg_shuff_pass),:)),...
%     'D:\Projects\PAE_PlaceCell\Figures\shufflepass\group1dg_shuff_pass','dpi','-r80')
% visualize_cells(uniqueRowsCA(group2dgid(logical(group2dg_shuff_pass),:)),...
%     'D:\Projects\PAE_PlaceCell\Figures\shufflepass\group2dg_shuff_pass','dpi','-r80')

%% interneuron filter

% [group1ca1_in,group1ca1id_in]=interneuron_filter(group1ca1,group1ca1id,varnames);
% [group2ca1_in,group2ca1id_in]=interneuron_filter(group2ca1,group2ca1id,varnames);
% [group1ca3_in,group1ca3id_in]=interneuron_filter(group1ca3,group1ca3id,varnames);
% [group2ca3_in,group2ca3id_in]=interneuron_filter(group2ca3,group2ca3id,varnames);
% [group1dg_in,group1dgid_in]=interneuron_filter(group1dg,group1dgid,varnames);
% [group2dg_in,group2dgid_in]=interneuron_filter(group2dg,group2dgid,varnames);
% 
% fig=figure;fig.Color=[1 1 1];
% AllStatsca1=stat_plot(group1ca1_in,group2ca1_in,{'Sacc','PAE'},varnames)
% fig=figure;fig.Color=[1 1 1];
% AllStatsca3=stat_plot(group1ca3_in,group2ca3_in,{'Sacc','PAE'},varnames)
% fig=figure;fig.Color=[1 1 1];
% AllStatsdg=stat_plot(group1dg_in,group2dg_in,{'Sacc','PAE'},varnames)

%% PLACE CELL FILTER
[group1ca1,group1ca1id]=placefieldfilter(group1ca1,group1ca1id,varnames);
[group2ca1,group2ca1id]=placefieldfilter(group2ca1,group2ca1id,varnames);
[group1ca3,group1ca3id]=placefieldfilter(group1ca3,group1ca3id,varnames);
[group2ca3,group2ca3id]=placefieldfilter(group2ca3,group2ca3id,varnames);
[group1dg,group1dgid]=placefieldfilter(group1dg,group1dgid,varnames);
[group2dg,group2dgid]=placefieldfilter(group2dg,group2dgid,varnames);



%% sessions for decoding
% decode=[unique(group1ca1id(:,1));...
% unique(group2ca1id(:,1));...
% unique(group1ca3id(:,1));...
% unique(group2ca3id(:,1));...
% unique(group1dgid(:,1));...
% unique(group2dgid(:,1))];
% save('D:\Projects\PAE_PlaceCell\decoding\decode','decode')


%%
[uCA,~,~] = uniqueRowsCA(group1ca1id);
disp([num2str(size(uCA,1)),' control ca1 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca1id);
disp([num2str(size(uCA,1)),' pae ca1 place cells'])

[uCA,~,~] = uniqueRowsCA(group1ca3id);
disp([num2str(size(uCA,1)),' control ca3 place cells'])
[uCA,~,~] = uniqueRowsCA(group2ca3id);
disp([num2str(size(uCA,1)),' pae ca3 place cells'])

[uCA,~,~] = uniqueRowsCA(group1dgid);
disp([num2str(size(uCA,1)),' control dg place cells'])
[uCA,~,~] = uniqueRowsCA(group2dgid);
disp([num2str(size(uCA,1)),' pae dg place cells'])

disp('...')
disp('control rats with ca1')
C=unique(extractBefore(unique(group1ca1id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group1ca1id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('pae rats with ca1')
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
disp('pae rats with ca3')
C=unique(extractBefore(unique(group2ca3id(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group2ca3id);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('control rats with dg')
C=unique(extractBefore(unique(group1dgid(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group1dgid);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end

disp('...')
disp('pae rats with dg')
C=unique(extractBefore(unique(group2dgid(:,1)),'_'));
for i=1:length(C)
    temp=uniqueRowsCA(group2dgid);
    disp([C{i},' ',num2str(sum(contains(temp(:,1),C{i})))])
end


%% compile and save as csv
id=[group1ca1id;group1ca3id;group1dgid;group2ca1id;group2ca3id;group2dgid];


id = [id,[cellstr(repmat('control',size([group1ca1id;group1ca3id;group1dgid],1),1));...
    cellstr(repmat('pae',size([group2ca1id;group2ca3id;group2dgid],1),1))]]

Rdata = [group1ca1;group1ca3;group1dg;group2ca1;group2ca3;group2dg];

Rdata = [num2cell(Rdata),id];

varnames = [varnames,{'session','tt','cell','area','group'}];
varnames = regexprep(varnames, '\W', '');
Rdata = cell2table(Rdata,'VariableNames',varnames);


writetable(Rdata,'F:\Projects\PAE_PlaceCell\Rdata_pae_sfn2019_lineartrack.csv')
%%

%  [pontential]=COM(group1ca1id)
% visualizecells(uniqueRowsCA(group1ca1id),'sacca1')
% visualizecells(uniqueRowsCA(group2ca1id),'paeca1')
% visualizecells(uniqueRowsCA(group1ca3id),'sacca3')
% visualizecells(uniqueRowsCA(group2ca3id),'paeca3')

fig=figure;fig.Color=[1 1 1];
AllStatsca1=stat_plot(group1ca1,group2ca1,{'Sacc','PAE'},varnames)
fig=figure;fig.Color=[1 1 1];
AllStatsca3=stat_plot(group1ca3,group2ca3,{'Sacc','PAE'},varnames)
fig=figure;fig.Color=[1 1 1];
AllStatsdg=stat_plot(group1dg,group2dg,{'Sacc','PAE'},varnames)



for i=1:length(varnames)
    fig=figure;fig.Color=[1 1 1];
    plot(sort(group1ca1(:,i)),linspace(0,1,length(group1ca1(:,i))),'LineWidth',2);hold on
    plot(sort(group1ca3(:,i)),linspace(0,1,length(group1ca3(:,i))),'LineWidth',2)
    plot(sort(group1dg(:,i)),linspace(0,1,length(group1dg(:,i))),'LineWidth',2)
    
    plot(sort(group2ca1(:,i)),linspace(0,1,length(group2ca1(:,i))),'LineWidth',2)
    plot(sort(group2ca3(:,i)),linspace(0,1,length(group2ca3(:,i))),'LineWidth',2)
    plot(sort(group2dg(:,i)),linspace(0,1,length(group2dg(:,i))),'LineWidth',2)
    
    xlabel(varnames{i})
    legend('control ca1','control ca3','control dg','pae ca1','pae ca3','pae dg','Location','best')
    grid on
    pause(.000001)
    export_fig(fullfile('F:\Projects\PAE_PlaceCell\Figures\region_comparison',[varnames{i},'.png']))

    close(fig)
end





ph1=group1ca3(group1ca3(:,contains(varnames,'Phpval'))<.05,contains(varnames,'PhcircLinCorr'));
ph2=group2ca3(group2ca3(:,contains(varnames,'Phpval'))<.05,contains(varnames,'PhcircLinCorr'));
figure
PHca3=stat_plot(ph1,ph2,{'Sacc','PAE'},varnames{contains(varnames,'PhcircLinCorr')},1)
fig=figure;fig.Color=[1 1 1]
[h1,h2]=CoolHistogram(ph1,ph2,50,varnames{contains(varnames,'PhcircLinCorr')})


for i=1:length(varnames)
    fig=figure('Name',['ca3 ',varnames{i}],'NumberTitle','off');
    AllStatsca3=stat_plot(group1ca3(:,i),group2ca3(:,i),{'Sacc','PAE'},varnames{i},'plots',2);
    toPPT(fig,'exportMode','matlab');
    toPPT('setTitle',AllStatsca3);
    close all
end

visualize_cells(uniqueRowsCA(group1ca1id),...
    'F:\Projects\PAE_PlaceCell\Figures\PlaceCells\control_ca1');
visualize_cells(uniqueRowsCA(group2ca1id),...
    'F:\Projects\PAE_PlaceCell\Figures\PlaceCells\pae_ca1');
visualize_cells(uniqueRowsCA(group1ca3id),...
    'F:\Projects\PAE_PlaceCell\Figures\PlaceCells\control_ca3');
visualize_cells(uniqueRowsCA(group2ca3id),...
    'F:\Projects\PAE_PlaceCell\Figures\PlaceCells\pae_ca3');
visualize_cells(uniqueRowsCA(group1dgid),...
    'F:\Projects\PAE_PlaceCell\Figures\PlaceCells\control_dg');
visualize_cells(uniqueRowsCA(group2dgid),...
    'F:\Projects\PAE_PlaceCell\Figures\PlaceCells\pae_dg');
% 
% visualizecells(group1ca1id,'control_ca1')
% visualizecells(group2ca1id,'pae_ca1')
% visualizecells(group1ca3id,'control_ca3')
% visualizecells(group2ca3id,'pae_ca3')


popvector(group1ca1id,group1ca1(:,end),'controlca1');
popvector(group1ca3id,group1ca3(:,end),'controlca3');

popvector(group2ca1id,group2ca1(:,end),'paeca1');
popvector(group2ca3id,group2ca3(:,end),'paeca3');


% AllStatsca1=stat_plot(group1ca1(:,contains(varnames,'thetaindex')),group2ca1(:,contains(varnames,'thetaindex')),{'Sacc','PAE'},'thetaindex',1)



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
%     AllStatsca3=stat_plot(group1ca3(randi([1,size(group1ca3,1)],13,1),:),group2ca3,{'Sacc','PAE'},varnames,3);
%     p=[ p,str2double(extractBetween(AllStatsca3,'p=',','))];
%     catch
%     end
% end
% (sum(p<.05,2)/size(p,2))>.95

AllStatsca1=stat_plot(group1ca1(:,1),group2ca1(:,1),{'Sacc','PAE'},varnames{1},2)
AllStatsca1=stat_plot(group1ca1(:,2),group2ca1(:,2),{'Sacc','PAE'},varnames{2},2)
AllStatsca1=stat_plot(group1ca1(:,4),group2ca1(:,4),{'Sacc','PAE'},varnames{4},2)
AllStatsca1=stat_plot(group1ca1(:,5),group2ca1(:,5),{'Sacc','PAE'},varnames{5},2)

AllStatsca1=stat_plot(group1ca1(:,19),group2ca1(:,19),{'Sacc','PAE'},varnames{19},2)

AllStatsca1=stat_plot(group1ca1(group1ca1(:,20)<.05,19),group2ca1(group2ca1(:,20)<.05,19),{'Sacc','PAE'},varnames{19},2)

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


stat_plot(nanstd(group1ca1_stab.com,0,2)/sqrt(size(group1ca1_stab.com,1)),...
    nanstd(group2ca1_stab.com,0,2)/sqrt(size(group2ca1_stab.com,1)),{'Sacc','PAE'},{'std_com'},2)

stat_plot(nanstd(group1ca3_stab.com,0,2)/sqrt(size(group1ca3_stab.com,1)),...
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

idx=group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'FieldWidth'))>=8 &...
    group(:,contains(varnames,'FieldWidth'))<=80 &...
    group(:,contains(varnames,'nlaps'))>=15 &...
    group(:,contains(varnames,'nSpikes'))>=100 &...
    group(:,contains(varnames,'InformationContent'))>=.15;

groupid=groupid(idx,:);

group=group(idx,:);
end

function [group,groupid]=interneuron_filter(group,groupid,varnames)
idx=group(:,contains(varnames,'spikewidth'))<=.25;
groupid=groupid(idx,:);

group=group(idx,:);
end

function [group,groupid]= preshuffle_placefieldfilter(group,groupid,varnames)
groupid=groupid(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'InformationContent'))>=0.25 &...
    group(:,contains(varnames,'nlaps'))>=10 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);

group=group(group(:,contains(varnames,'PeakRate'))>=1 &...
    group(:,contains(varnames,'InformationContent'))>=0.25 &...
    group(:,contains(varnames,'nlaps'))>=10 &...
    group(:,contains(varnames,'nSpikes'))>=100,:);
end

function visualizecells(groupid,group)
cd('D:\Projects\PAE_PlaceCell\ProcessedData')
for i=1:length(groupid)
    data=load(groupid{i,1});
    try
        postprocessFigures.main(data,{groupid{i,2},str2double(groupid(i,3))},'colorcode','HD');
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
    ['F:\Projects\PAE_PlaceCell\Figures\popvectors\',group,filesep,'popvec.png'])

% savefig(['F:\Projects\PAE_PlaceCell\Figures\popvectors\',group,filesep,'popvec.fig'])
close all
end