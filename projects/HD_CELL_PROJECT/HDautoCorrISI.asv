% HDautoCorrISI
% opens .ts_R text files containing HD cell data and plots autoCorrISI
%
%
% Ryan E Harvey 2018
%
% cd to data and get file names
clear;clc;close all
com=which('HDautoCorrISI');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'BClarkToolbox',filesep, 'Analysis']);

if ismac
    path='/Users/ryanharvey/Downloads/HeadDIrectionCells_LauraRyan';
else
    path='D:\Projects\Multi_Region_HD\HeadDIrectionCells_wRfiles';
end

cd(path)
% files=dir( '**/*.ts_R');
files=dir('**/*_Timestamps*.*');
for i=1:length(files)
    folders{i,:}=[files(i).folder,filesep,files(i).name];
end
disp(folders)

for i=1:length(folders)
    area=strsplit(folders{i},filesep);
    area=area{end};
    area=erase(area,'_Timestamps');
    areas{i}=erase(area,'_');
end
disp(areas)
clear area files 

for a=1:length(folders)
    cd(folders{a})
    filenames=dir;
    filenames={filenames.name}';
    filenames(strcmp(filenames,'..') | strcmp(filenames,'.') )=[];
    disp([num2str(length(filenames)),' ',areas{a},' cells'])
    for i=1:length(filenames)
        % OPEN TEXT FILE AND CLEAN DATA
        disp(['Running:  ',filenames{i}])
        
        fileID=fopen(filenames{i},'r');
        dataArray = textscan(fileID, '%f%*s%[^\n\r]', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,2-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        ts=[dataArray{1:end-1}];
        clearvars fileID dataArray ans;
        
        ts(ts<0)=[];
        
        normISIsmooth(i,:)=autoCorrISI(ts,0);
        [thetaindex(i,1),peak(i,:),cor(i,:),lag] = thetamodulation((ts./100)./1000);
        
        
        %% Bin in 2ms bins
        spk=(ts./100)./1000;
        max_lag = 0.5;
        t_bin=0.002; % 2ms
        % Acor - taken from intrinsic frequency 2
        if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
            max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
        end
        [two_ms_cor, lag] = CrossCorr(spk, spk, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', 'prob');
        if length(lag)~=length(two_ms_cor)
            lag=linspace(-.5,.5,length(two_ms_cor));
        end
        
        two_ms_cor = two_ms_cor/max(two_ms_cor(lag~=0))*2-1;
        if any(two_ms_cor>1)
            two_ms_cor=rescale(two_ms_cor,-1,1);
        end
        two_ms_cor_save(i,:)=two_ms_cor;
        
        
        %% theta skip
        max_lag=.6;
        t_bin=0.005; % 5ms
        if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
            max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
        end
        [skipcor, lag] = CrossCorr(spk, .6, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', 'prob');
        thetaskip(i,1)=ThetaSkipping(skipcor,lag,t_bin);
        
        %% Bin in .5ms bins
        spk=(ts./100)./1000;
        max_lag = 0.5;
        t_bin=0.0005;% .5ms bins
        % Acor - taken from intrinsic frequency 2
        if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
            max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
        end
        [two_ms_cor, lag] = CrossCorr(spk, spk, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', 'prob');
        if length(lag)~=length(two_ms_cor)
            lag=linspace(-.5,.5,length(two_ms_cor));
        end
        two_ms_cor = two_ms_cor/max(two_ms_cor(lag~=0))*2-1;
        if any(two_ms_cor>1)
            two_ms_cor=rescale(two_ms_cor,-1,1);
        end
        half_ms_cor_save(i,:)=two_ms_cor;
        
        
        
        % Burst Index: (Coletta et al. 2018)
         [bi(i,:),isi(i,:),isilag]=burstindex((ts./100)./1000,.0005);
        % 
        
%         x=linspace(ts(1),ts(end),length(ts));
%         thetashuff=zeros(500,1);
%         parfor ishuff=1:500
%             [thetashuff(ishuff,:),~,~,~]=thetamodulation((sort(datasample(x,length(ts)))./100)./1000);
%         end
%         modulated(i,1)=thetaindex(i)>=prctile(thetashuff,95);
%         
        ids{i,1}=filenames{i};
    end
    HDdata.(areas{a}).id=ids;
    HDdata.(areas{a}).thetaindex=thetaindex;
    HDdata.(areas{a}).peak=peak;
    HDdata.(areas{a}).autocor=cor;
    HDdata.(areas{a}).two_ms_cor=two_ms_cor_save;
    HDdata.(areas{a}).half_ms_cor=half_ms_cor_save;
    HDdata.(areas{a}).BurstIndex=bi;
    HDdata.(areas{a}).isi=isi;
    HDdata.(areas{a}).thetaskip=thetaskip;

    

%     HDdata.(areas{a}).thetamod=modulated;

    clear thetaindex peak cor normISIsmooth modulated ids two_ms_cor_save half_ms_cor_save bi isi thetaskip
    
end

%% PLOT

% theta skipping by layer
clear data
data.ATN=HDdata.ATN.thetaskip;
data.PoSDeep=HDdata.PoSDeep.thetaskip;
data.PoSSup=HDdata.PoSSup.thetaskip;
data.MECDeepLayers=HDdata.MECDeepLayers.thetaskip;
data.MECSupLayers=HDdata.MECSupLayers.thetaskip;
data.PaSDeep=HDdata.PaSDeep.thetaskip;
data.PaSSup=HDdata.PaSSup.thetaskip;

fig=figure;fig.Color=[1 1 1];
ECDF_plot(data,'Theta Skipping');



% close all
mec=[HDdata.MECSupLayers.thetamod;HDdata.MECDeepLayers.thetamod];
pas=[HDdata.PaSSup.thetamod;HDdata.PaSDeep.thetamod];
pos=[HDdata.PoSSup.thetamod;HDdata.PoSDeep.thetamod];
atn=[HDdata.ATN.thetamod];

y=[sum(mec==1)/length(mec),sum(mec==0)/length(mec);...
    sum(pas==1)/length(pas),sum(pas==0)/length(pas);...
    sum(pos==1)/length(pos),sum(pos==0)/length(pos);...
    sum(atn==1)/length(atn),sum(atn==0)/length(atn)];
categories={'MEC','PaS','PoS','ATN'};

barfigregion=figure;barfigregion.Color=[1 1 1];
b=barh(y,'stacked')
xlabel('Proportion')
ylabel('Layer/Region')
set(b,'LineStyle','none');
set(gca,'yticklabel',categories,'box','off','LineWidth',2,'FontWeight','bold','FontSize',20,'TickLength',[0;0])
b(1, 1).FaceColor   = [0.929,  0.694,  0.125];
b(1, 2).FaceColor   = [0.3,  0.3,  0.3];
%%

% theta index by layer
clear data
data.ATN=HDdata.ATN.thetaindex;
data.PoSDeep=HDdata.PoSDeep.thetaindex;
data.PoSSup=HDdata.PoSSup.thetaindex;
data.MECDeepLayers=HDdata.MECDeepLayers.thetaindex;
data.MECSupLayers=HDdata.MECSupLayers.thetaindex;
data.PaSDeep=HDdata.PaSDeep.thetaindex;
data.PaSSup=HDdata.PaSSup.thetaindex;

fig=figure;fig.Color=[1 1 1];
ECDF_plot(data,'Theta Index');
%% theta index by region
clear data
data.ATN=HDdata.ATN.thetaindex;
data.PoS=[HDdata.PoSDeep.thetaindex;HDdata.PoSSup.thetaindex];
data.MEC=[HDdata.MECDeepLayers.thetaindex;HDdata.MECSupLayers.thetaindex];
data.PaS=[HDdata.PaSDeep.thetaindex;HDdata.PaSSup.thetaindex];

fig=figure;fig.Color=[1 1 1];
ECDF_plot(data,'Theta Index');


%% ATN pop vec
[thetaindex_sorted,I]=sort(HDdata.ATN.thetaindex);
corSorted=HDdata.ATN.autocor(I,:);

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(corSorted);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('ATN')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(thetaindex_sorted(round(linspace(1,size(corSorted,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,201,2),'XTickLabel',...
    [-500 500],'FontSize',20,'FontWeight','bold','TickLength',[0,0])
%% PoS pop vec
[thetaindex_sorted,I]=sort([HDdata.PoSDeep.thetaindex;HDdata.PoSSup.thetaindex]);
tempcor=[HDdata.PoSDeep.autocor;HDdata.PoSSup.autocor];
corSorted=tempcor(I,:);

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(corSorted);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('PoS')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(thetaindex_sorted(round(linspace(1,size(corSorted,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,201,2),'XTickLabel',...
    [-500 500],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

%% MEC pop vec
[thetaindex_sorted,I]=sort([HDdata.MECDeepLayers.thetaindex;HDdata.MECSupLayers.thetaindex]);
tempcor=[HDdata.MECDeepLayers.autocor;HDdata.MECSupLayers.autocor];
corSorted=tempcor(I,:);

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(corSorted);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('MEC')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(thetaindex_sorted(round(linspace(1,size(corSorted,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,201,2),'XTickLabel',...
    [-500 500],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

%% PaS pop vec
[thetaindex_sorted,I]=sort([HDdata.PaSDeep.thetaindex;HDdata.PaSSup.thetaindex]);
tempcor=[HDdata.PaSDeep.autocor;HDdata.PaSSup.autocor];
corSorted=tempcor(I,:);

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(corSorted);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('PaS')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(thetaindex_sorted(round(linspace(1,size(corSorted,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,201,2),'XTickLabel',...
    [-500 500],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

%% ////////////////////// 60ms AUTO CORRS ///////////////////////////////////
%% ATN pop vec
[thetaindex_sorted,I]=sort(HDdata.ATN.thetaindex);
corSorted=HDdata.ATN.autocor(I,51:57);

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(corSorted);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('ATN')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(thetaindex_sorted(round(linspace(1,size(corSorted,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,7,2),'XTickLabel',...
    [0 60],'FontSize',20,'FontWeight','bold','TickLength',[0,0])
%% PoS pop vec
[thetaindex_sorted,I]=sort([HDdata.PoSDeep.thetaindex;HDdata.PoSSup.thetaindex]);
tempcor=[HDdata.PoSDeep.autocor;HDdata.PoSSup.autocor];
corSorted=tempcor(I,51:57);

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(corSorted);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('PoS')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(thetaindex_sorted(round(linspace(1,size(corSorted,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,7,2),'XTickLabel',...
    [0 60],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

%% MEC pop vec
[thetaindex_sorted,I]=sort([HDdata.MECDeepLayers.thetaindex;HDdata.MECSupLayers.thetaindex]);
tempcor=[HDdata.MECDeepLayers.autocor;HDdata.MECSupLayers.autocor];
corSorted=tempcor(I,51:57);

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(corSorted);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('MEC')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(thetaindex_sorted(round(linspace(1,size(corSorted,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,7,2),'XTickLabel',...
    [0 60],'FontSize',20,'FontWeight','bold','TickLength',[0,0])
%% PaS pop vec
[thetaindex_sorted,I]=sort([HDdata.PaSDeep.thetaindex;HDdata.PaSSup.thetaindex]);
tempcor=[HDdata.PaSDeep.autocor;HDdata.PaSSup.autocor];
corSorted=tempcor(I,51:57);

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(corSorted);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('PaS')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(thetaindex_sorted(round(linspace(1,size(corSorted,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,7,2),'XTickLabel',...
    [0 60],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

%% \\\\\\\\\\\\\\\\\\\\\\\\\\60sec autocor using 2ms bins \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%%

for i=1:length(areas)
tempcor=[HDdata.(areas{i}).two_ms_cor(:,251:281)];

[s,I]=sort(HDdata.(areas{i}).BurstIndex,'descend');

tempcor=tempcor(I,:);
% [tempcor]=arrangenorm(tempcor);


fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[-6,34,476,924];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(tempcor);
box off;
ylabel('Cells')
xlabel('Lag (ms)')
title(areas{i})

set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,size(tempcor,2),2),'XTickLabel',...
    [0 60],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

end


idx(1)=find(lag==0);
idx2=find(lag>=0.012);
idx(2)=idx2(1);

forpca=[];
for i=1:length(areas)
    forpca=[forpca;HDdata.(areas{i}).half_ms_cor(:,idx(1):idx(2))];
end

% forpca=[HDdata.MECSupLayers.half_ms_cor(:,idx(1):idx(2));HDdata.MECDeepLayers.half_ms_cor(:,idx(1):idx(2))];

[forpca]=arrangenorm(forpca);
[coeff,score,latent,tsquared,explained]=pca(forpca);
figure;scatter3(score(:,1),score(:,2),score(:,3),'Filled','k');hold on
xlabel('PCA 1')
ylabel('PCA 2')
zlabel('PCA 3')

[idx,C]=kmeans(score,2);

scatter3(score(idx==1,1),score(idx==1,2),score(idx==1,3),'Filled','k');hold on
scatter3(score(idx==2,1),score(idx==2,2),score(idx==2,3),'Filled','r');hold on

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[-6,34,476,924];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(forpca(idx==1,:));
box off;
ylabel('Cells')
xlabel('Lag (ms)')
title('kmeans 1')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,size(forpca,2),2),'XTickLabel',...
    [0 12],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[-6,34,476,924];
set(fig,'defaultAxesColorOrder',[0 0 0]);
imagesc(forpca(idx==2,:));
box off;
ylabel('Cells')
xlabel('Lag (ms)')
title('kmeans 2')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,size(forpca,2),2),'XTickLabel',...
    [0 12],'FontSize',20,'FontWeight','bold','TickLength',[0,0])


%% PLOT and SAVE ISI

idx(1)=find(isilag==0);
idx2=find(isilag>=0.012);
idx(2)=idx2(1);

for i=1:length(areas)
    for ii=1:size(HDdata.(areas{i}).isi,1)
        tempcor=HDdata.(areas{i}).isi(ii,idx(1):idx(2));
        bi=HDdata.(areas{i}).BurstIndex(ii,1);
        x=linspace(0,12,length(tempcor));
        fig=figure;fig.Color=[1 1 1];
        b=bar(x,tempcor,'histc');
        b.FaceColor='k';
        b.EdgeColor='k';
        box off;
        axis tight
        ylabel('Spk Count')
        xlabel('Lag (ms)')
        id=strsplit(HDdata.(areas{i}).id{ii},'.');
        id=id{1};
        disp(id)
        title([id,' Burst Index: ',num2str( bi)])
        set(gca,'FontSize',10,'FontWeight','bold','TickLength',[0,0])
        
        cd(path)
        cd ..
        if ~exist([cd,filesep,areas{i},'_isi_Figures'],'dir')
            mkdir([cd,filesep,areas{i},'_isi_Figures'])
        end
        print(fig,'-dpng', '-r150',[cd,filesep,areas{i},'_isi_Figures',filesep,id,'.png'])
        close all
    end
end

%% Plot and Save Autocorr


for i=1:length(areas)
    for ii=1:size(HDdata.(areas{i}).autocor,1)
        tempcor=HDdata.(areas{i}).autocor(ii,:);
        
        bi=HDdata.(areas{i}).BurstIndex(ii,1);
        thetaindex=HDdata.(areas{i}).thetaindex(ii,1);

        fig=figure;fig.Color=[1 1 1];
        plot(tempcor,'LineWidth',2, 'color','k');
        axis tight
        hold on;box off; axis off

        id=strsplit(HDdata.(areas{i}).id{ii},'.');
        id=id{1};
        disp(id)
        
        title([' BurstIndex: ',num2str(bi),' ThetaIndex: ',num2str(thetaindex)])
 
        if ~exist(['D:\Projects\Multi_Region_HD\autocor',filesep,areas{i}],'dir')
            mkdir(['D:\Projects\Multi_Region_HD\autocor',filesep,areas{i}])
        end
        print(fig,'-dpng', '-r150',['D:\Projects\Multi_Region_HD\autocor',filesep,areas{i},filesep,id,'.png'])
        close all
    end
end

