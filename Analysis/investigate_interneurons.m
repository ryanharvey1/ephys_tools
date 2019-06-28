% investigate_interneurons

% Re: interneuron vs. excitatory cells - could you make a plot of burst index vs. peak to trough as in Fig 1A of the Buszaki paper?
% 
% Also, could you make a separate plot for each region and a plot pooled across the three hippocampal regions?

control={'RH13','RH14','LS21','LS23','LE2821','LE2823','LEM3116','LEM3120'};
pae={'RH11','RH16','LS17','LS19','LE2813','LEM3124'};

cd D:\Projects\PAE_PlaceCell\ProcessedData
sessions=dir('*.mat');

c=1;
for i=1:length(sessions)
    load(fullfile(sessions(i).folder,sessions(i).name),'Spikes','spikesID','avgwave','frames');
    disp(sessions(i).name)
    for s=1:length(Spikes)
        disp(['cell',num2str(c)])
        
        avgrate(c,1)=length(Spikes{s})/frames(end,1);
        
        [burstIdx(c,1)]=burst_index(Spikes{s});
        
        prop=waveformprop(avgwave(s));
        spikewidth(c,1)=prop(3);
        
        [p,I]=max(max(avgwave{s},[],2));
        peakwave(c,:)=rescale(avgwave{s}(I,:),0,1);
        
        id(c,:)=[sessions(i).name,spikesID.TetrodeNum(s),num2str(spikesID.CellNum(s))];
        c=c+1;
    end
end

%%

% snaped_sessions=snapped;
% for i=1:length(snaped_sessions)
%     ratID=strsplit(snaped_sessions{i},filesep);
%     rat=ratID{end-1};
%     sessionID=['S',strjoin(regexp(ratID{end},'\d*','Match'),'')];
%     snaped_sessions{i}=[rat,'_',sessionID];
% end
% snap_idx=contains({sessions.name}',snaped_sessions);
% 
% all_peakwave=peakwave;
% for ii=0:1
% peakwave=all_peakwave(snap_idx==ii,:);
% 
% figure(ii+1)
% subplot(3,3,2.5)
% [coeff,score,latent,tsquared,explained,mu] = pca(peakwave);
% scatter3(score(:,1),score(:,2),score(:,3),'k','Filled')
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')
% title('First 3 Components')
% 
% clust = zeros(size(score,1),10);
% for i=1:10
%     clust(:,i) = kmeans(score,i,'emptyaction','singleton',...
%         'replicate',5);
% end
% eva = evalclusters(score,clust,'CalinskiHarabasz');
% 
% 
% subplot(3,3,4)
% plot(eva.CriterionValues,'k','LineWidth',2)
% hold on
% plot(eva.OptimalK,eva.CriterionValues(eva.OptimalK),'*r')
% ylabel('Criterion Values')
% xlabel('Number of Clusters')
% title('K-means cluster evaluation')
% box off
% grid on
% 
% subplot(3,3,5)
% Z = linkage(peakwave,'ward','euclidean');
% color = Z(end-eva.OptimalK+2,3)-eps;
% H = dendrogram(Z,0,'ColorThreshold',color);
% set(gca,'xtick',[])
% ylabel('Distance')
% title(['Dendrogram parsed at ',num2str(eva.OptimalK),' clusters'])
% 
% % recolor lines to match cluster color
% % lineColours = cell2mat(get(H,'Color'));
% % colourList = unique(lineColours, 'rows');
% % myColours = [0 0 1;1 0 0];
% % for colour = 2:size(colourList,1)
% %     idx = ismember(lineColours, colourList(colour,:), 'rows');
% %     lineColours(idx, :) = repmat(myColours(colour-1,:),sum(idx),1);
% % end
% % for line = 1:size(H,1)
% %     set(H(line), 'Color', lineColours(line,:));
% % end
% 
% T = cluster(Z,'maxclust',eva.OptimalK);
% 
% 
% subplot(3,3,6)
% scatter3(score(:,1),score(:,2),score(:,3),20,T,'Filled');
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')
% title('First 3 Components Clustered')
% colormap(gca,jet) 
% 
% 
% subplot(3,3,7)
% imagesc(peakwave(T==1,:));
% title('Cluster 1 (blue cluster)')
% box off
% colormap(gca,viridis(255))
% xlabel('Microseconds')
% ylabel('Waveforms')
% set(gca,'XTick',linspace(1,size(peakwave,2),4),'XTicklabel',sprintf('   %.1f\n',linspace(1,size(peakwave,2),4)/100))
% 
% subplot(3,3,8)
% imagesc(peakwave(T==2,:))
% title('Cluster 2 (red cluster)')
% box off
% colormap(gca,viridis(255))
% xlabel('Microseconds')
% ylabel('Waveforms')
% set(gca,'XTick',linspace(1,size(peakwave,2),4),'XTicklabel',sprintf('   %.1f\n',linspace(1,size(peakwave,2),4)/100))
% 
% 
% subplot(3,3,9)
% Group1=peakwave(T==1,:);
% Group2=peakwave(T==2,:);
% SEM=nanstd(Group1)/sqrt(size(Group1,1));
% h1=shadedErrorBar(1:size(Group1,2),nanmean(Group1),SEM,'-b',0);
% hold on
% SEM=nanstd(Group2)/sqrt(size(Group2,1));
% h2=shadedErrorBar(1:size(Group2,2),nanmean(Group2),SEM,'-r',0);
% xlim([1 size(Group1,2)])
% set(gca,'Box','off','LineWidth',2,...
%     'XTick',linspace(1,size(peakwave,2),4),'XTicklabel',sprintf('   %.1f\n',linspace(1,size(peakwave,2),4)/100))
% xlabel('Microseconds')
% ylabel('Normalized amplitude')
% title('Average Waveforms')
% 
% end
%%
% get ids
id=get_region_id(id,'D:\Projects\PAE_PlaceCell\AnimalMetadata');

% get interneurons
interneuron_idx=spikewidth<.25;



burstIdxca1 = burstIdx(strcmp(id(:,4),'ca1'),:);
burstIdxca3 = burstIdx(strcmp(id(:,4),'ca3'),:);
burstIdxdg = burstIdx(strcmp(id(:,4),'dg'),:);

spikewidthca1 = spikewidth(strcmp(id(:,4),'ca1'),:);
spikewidthca3 = spikewidth(strcmp(id(:,4),'ca3'),:);
spikewidthdg = spikewidth(strcmp(id(:,4),'dg'),:);

avgrateca1 = avgrate(strcmp(id(:,4),'ca1'),:);
avgrateca3 = avgrate(strcmp(id(:,4),'ca3'),:);
avgratedg = avgrate(strcmp(id(:,4),'dg'),:);




fig=figure;fig.Color=[1 1 1];
subplot(1,4,1)
plot([spikewidthca1;spikewidthca3;spikewidthdg],[burstIdxca1;burstIdxca3;burstIdxdg],'.k')
set(gca, 'YScale', 'log')
ylabel('burst index')
xlabel('spike width')
title('all')
subplot(1,4,2)
plot(spikewidthca1,burstIdxca1,'.k')
set(gca, 'YScale', 'log')

ylabel('burst index')
xlabel('spike width')
title('ca1')
subplot(1,4,3)
plot(spikewidthca3,burstIdxca3,'.k')
set(gca, 'YScale', 'log')

ylabel('burst index')
xlabel('spike width')
title('ca3')
subplot(1,4,4)
plot(spikewidthdg,burstIdxdg,'.k')
set(gca, 'YScale', 'log')
ylabel('burst index')
xlabel('spike width')
title('dg')



%%
fig=figure;fig.Color=[1 1 1];
scatterhist([spikewidthca1;spikewidthca3;spikewidthdg],...
    log([burstIdxca1;burstIdxca3;burstIdxdg]),...
    'NBins',[100,100],'Kernel','overlay','Color',...
    'k','Marker','.','MarkerSize',20)
ylabel('burst index')
xlabel('spike width')
title('all')

fig=figure;fig.Color=[1 1 1];
subplot(1,4,2)
scatterhist(spikewidthca1,log(burstIdxca1),...
    'NBins',[100,100],'Kernel','overlay','Color',...
    'k','Marker','.','MarkerSize',20)
ylabel('burst index')
xlabel('spike width')
title('ca1')

fig=figure;fig.Color=[1 1 1];
subplot(1,4,3)
scatterhist(spikewidthca3,log(burstIdxca3),...
    'NBins',[100,100],'Kernel','overlay','Color',...
    'k','Marker','.','MarkerSize',20)
ylabel('burst index')
xlabel('spike width')
title('ca3')

fig=figure;fig.Color=[1 1 1];
subplot(1,4,4)
scatterhist(spikewidthdg,log(burstIdxdg),...
    'NBins',[100,100],'Kernel','overlay','Color',...
    'k','Marker','.','MarkerSize',20)
ylabel('burst index')
xlabel('spike width')
title('dg')

%%
fig=figure;fig.Color=[1 1 1];
subplot(1,4,1)
scatterhist([spikewidthca1;spikewidthca3;spikewidthdg],...
    [avgrateca1;avgrateca3;avgratedg],...
    'NBins',[100,100],'Kernel','on','Color',...
    'k','Marker','.','MarkerSize',20)
ylabel('avg rate (hz)')
xlabel('spike width')
title('all')

fig=figure;fig.Color=[1 1 1];
subplot(1,4,2)
scatterhist(spikewidthca1,avgrateca1,...
    'NBins',[100,100],'Kernel','on','Color',...
    'k','Marker','.','MarkerSize',20)
ylabel('avg rate (hz)')
xlabel('spike width')
title('ca1')

fig=figure;fig.Color=[1 1 1];

subplot(1,4,3)
scatterhist(spikewidthca3,avgrateca3,...
    'NBins',[100,100],'Kernel','on','Color',...
    'k','Marker','.','MarkerSize',20)
ylabel('avg rate (hz)')
xlabel('spike width')
title('ca3')

fig=figure;fig.Color=[1 1 1];

subplot(1,4,4)
scatterhist(spikewidthdg,avgratedg,...
    'NBins',[100,100],'Kernel','on','Color',...
    'k','Marker','.','MarkerSize',20)
ylabel('avg rate (hz)')
xlabel('spike width')
title('dg')

function [bi]=burst_index(spk,t_bin)
% burst index = 3-5msec/200-300msec

% Input: spk: timestamps in seconds
% Output: bi: Burst Index

% Ryan Harvey

max_lag = 0.5;
if ~exist('t_bin','var')
    t_bin=0.001; 
end
% Acor - taken from intrinsic frequency 2
if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
    max_lag = t_bin*floor(max_lag/t_bin)+max_lag*t_bin;
end
% tic
% [cor, lag] = CrossCorr(spk, spk, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', 'count');
% toc
% bi=sum(cor(504:506))/sum(cor(701:801));




[cor,lag] = MClustStats.AutoCorr(spk, .001, 501);

bi=nansum(cor(lag>0.003 & lag<0.005))/nansum(cor(lag>0.2 & lag<0.3));
if isnan(bi)
    bi=0;
end
end

function snaped_sessions=snapped()
sessions = dir('D:\Projects\PAE_PlaceCell\data\**\*Sorted');
for s=1:length(sessions)
    tts=dir(fullfile(sessions(s).folder,sessions(s).name,'*info.mat'));
    load(fullfile(tts(1).folder,tts(1).name));
    snap(s,1)=~any(isnan(confidence)) & ~any(isnan(final_grades));
end
snaped_sessions={sessions(snap).folder}';

end
