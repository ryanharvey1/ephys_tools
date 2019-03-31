% AnalyzeResults_ATNad

addpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\Visualize')
data=compileResults('F:\ClarkP30_Recordings\ProcessedData');

control={'ATN04','ATN05','ATN09','ATN08','LB03'};%,'LB01','LB03','LB05'};
transgenic={'ATN01','ATN03','ATN07','LB04'}; %,'LB02','LB04','LB06'};

%% COMPILE GROUPS
data.control.measures=[];
data.control.id=[];
for i=1:length(control)
    data.control.measures=cat(1,data.control.measures,data.(control{i}).measures);
    data.control.id=cat(1,data.control.id,data.(control{i}).id);
end

data.transgenic.measures=[];
data.transgenic.id=[];
for i=1:length(transgenic)
    data.transgenic.measures=cat(1,data.transgenic.measures,data.(transgenic{i}).measures);
    data.transgenic.id=cat(1,data.transgenic.id,data.(transgenic{i}).id);
end


% COMPILE DATA
group1=data.control.measures;
group2=data.transgenic.measures;
group1id=data.control.id;
group2id=data.transgenic.id;


%% DELETE MEASURES FOR LINEAR TRACK
varnames=data.varnames;
group1(isinf(group1))=NaN;
group2(isinf(group2))=NaN;

% thetaFreq=CDFplots(group1(:,23,1),group2(:,23,1),{'Control','TgF344-AD'},varnames{23},1)
% thetaFreq=CDFplots(group1(group1(:,24,1)<2500000,24,1),group2(:,24,1),{'Control','TgF344-AD'},varnames{24},1

colstodelete=contains(varnames,["Displacement", "DisplacementCorr","DirectionalityIndex","nlaps","rateoverlap"...
    ,"fieldoverlap","lap_perm_stability","stabilityoverlaps","meanstability","spatialcorrelation"]);

varnames(colstodelete)=[];
group1(:,colstodelete,:)=[];
group2(:,colstodelete,:)=[];


%% PLACE CELL FILTER
[group1place,group1placeid]=placefieldfilter(group1(:,:,:),group1id(:,:,1),varnames);
[group2place,group2placeid]=placefieldfilter(group2(:,:,:),group2id(:,:,1),varnames);

for i=1:length(varnames)
    
    all=CDFplots(group1place(:,i,1),group2place(:,i,1),{'Control','TgF344-AD'},varnames{i},2)
    toPPT(gcf,'exportMode','matlab'); % default would be "exportFig" instead of "matlab"
    toPPT('setTitle',all);
    
    close all
    
end

% visualizecells(group1placeid,'WTplace')
% visualizecells(group2placeid,'TGplace')

% [h1,h2]=CoolHistogram(group1place(:,1),group2place(:,1),100,'InfoContent')
% 
group1place(:,24,3)=abs(wrapTo180(group1place(:,24,3)));
group2place(:,24,3)=abs(wrapTo180(group2place(:,24,3)));

AllStats=CDFplots(group1place(:,:,3),group2place(:,:,3),{'Control','TgF344-AD'},varnames,2)

% done=PlotsStat_ATN(group1place(:,:,1:4),group2place(:,:,1:4),varnames)

g1tuning=get_tuning(group1placeid);
[g1tuning]=arrangenorm(g1tuning);
fields=fieldnames(g1tuning);
h = fspecial('gaussian');
fig=figure;
 fig.Color=[1 1 1];
    p = panel(fig); 
    p.pack(2, 4); 
    p.de.margin = 4;
   
for i=1:length(fields)
p(1, i).select(); 
ax1=gca; 
imagesc(filter2(h, g1tuning.(fields{i}))); colormap(ax1,parula);
axis off;axis tight
end
for i=1:4
p(2, i).select(); 
ax2=gca; 
imagesc(g1tuning.(fields{i})==1); colormap(ax2,flipud(gray));
axis off;axis tight
end


g2tuning=get_tuning(group2placeid);
[g2tuning]=arrangenorm(g2tuning);
fields=fieldnames(g2tuning);
h = fspecial('gaussian');
fig=figure;
 fig.Color=[1 1 1];
    p = panel(fig); 
    p.pack(2, 4); 
    p.de.margin = 4;
for i=1:length(fields)
p(1, i).select(); 
ax1=gca; 
imagesc(filter2(h, g2tuning.(fields{i}))); colormap(ax1,parula);
axis off;axis tight
end
for i=1:4
p(2, i).select(); 
ax2=gca; 
imagesc(g2tuning.(fields{i})==1); colormap(ax2,flipud(gray));
axis off;axis tight
end


%% HD Filter 
[group1HD,group1HDid]=headdirectionfilter(group1(:,:,:),group1id(:,:,1),varnames);
[group2HD,group2HDid]=headdirectionfilter(group2(:,:,:),group2id(:,:,1),varnames);
% 
% visualizecells(group1HDid,'WTHD')
% visualizecells(group2HDid,'TGHD')

% [h1,h2]=CoolHistogram(group1HD(:,9),group2HD(:,9),100,'Mean Vector Length')
% 
% AllStats=CDFplots(group1HD(:,:,1),group2HD(:,:,1),{'Control','TgF344-AD'},varnames,1)



%% INTERNEURON FILTER
[group1inter,group1interid]=interneuronfilter(group1,group1id,varnames)
[group2inter,group2interid]=interneuronfilter(group2,group2id,varnames)

[h1,h2]=CoolHistogram(group1inter(:,34),group2inter(:,34),100,'thetaindex')

% AllStats=CDFplots(group1inter(:,:,1),group2inter(:,:,1),{'Control','TgF344-AD'},varnames,1)

group1autocorr=popVecAutocorr(group1interid);
[g1thetaindex_sorted,I]=sort(group1inter(:,34,1));
g1autocorrPopvec=group1autocorr(I,:);
g1autocorrPopvec(sum(isnan(g1autocorrPopvec),2)>0,:)=[];

group2autocorr=popVecAutocorr(group2interid);
[g2thetaindex_sorted,I]=sort(group2inter(:,34,1));
g2autocorrPopvec=group2autocorr(I,:);
g2autocorrPopvec(sum(isnan(g2autocorrPopvec),2)>0,:)=[];




% figure; subplot(1,2,1); imagesc(g1autocorrPopvec); 
% 
% subplot(1,2,2);   imagesc(g2autocorrPopvec); 

fig=figure;fig.Color=[1 1 1];
fig.OuterPosition=[1 6 960 1052];
set(fig,'defaultAxesColorOrder',[0 0 0]);
subplot(1,2,1);
imagesc(g1autocorrPopvec);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('WT')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(g1thetaindex_sorted(round(linspace(1,size(g1autocorrPopvec,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,101,2),'XTickLabel',...
    [-500 500],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

subplot(1,2,2);
imagesc(g2autocorrPopvec);
box off;axis xy
ylabel('Cells')
xlabel('Lag (ms)')
title('TgF344-AD')
yyaxis right
set(gca,'YTick',linspace(0,1,8),'YTickLabel',...
    round(g2thetaindex_sorted(round(linspace(1,size(g2autocorrPopvec,1),8)),1)',2))
ylabel('Theta Index')
set(gca,'XMinorTick','on','YMinorTick','off',...
    'LineWidth',1,'box','off','XTick',linspace(1,101,2),'XTickLabel',...
    [-500 500],'FontSize',20,'FontWeight','bold','TickLength',[0,0])

%% ___________________________LOCAL FUNCTION BELOW_________________________
function [group,groupid]=placefieldfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 2 Hz, 
% 2) Minimum field width of 8 cm, 
% 3) Maximum field width of 180 cm, 
% 4) Mean vector length less than .5 to limit direction x place. 
% 5) at least 100 spikes
%IDX HD cells
% idx=sum(group(:,contains(varnames,'PeakRate'),1)>=1,3)>0 &...
%     sum(group(:,contains(varnames,'mean_vector_length'),1)>=.2,3)>0 &...
%     sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0;

%IDX place cells
idx=sum(group(:,contains(varnames,'PeakRate'),1)>=1,3)>0 &...
    sum(group(:,contains(varnames,'FieldWidth'),1)>=8,3)>0 &...
    sum(group(:,contains(varnames,'FieldWidth'),1)<=40,3)>0 & ...
    sum(group(:,contains(varnames,'InformationContent'),1)>=.4,3)>0 & ...
    sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);

end

%%
function [group,groupid]=headdirectionfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 2 Hz, 
% 2) Maximum field width of 8 cm,  
% 3) at least 100 spikes
idx=sum(group(:,contains(varnames,'PeakRate'),1)>=2,3)>0 &...
    sum(group(:,contains(varnames,'mean_vector_length'),1)>=.5,3)>0 &...
    sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);
end


function [group,groupid]=qualityfilter(group,groupid,varnames)
% 1) Minimum peak firing rate of 2 Hz, 
% 2) Maximum field width of 8 cm,  
% 3) at least 100 spikes
idx=sum(group(:,contains(varnames,'ShortISI'),1)<.05,3)>0 &...
    sum(group(:,contains(varnames,'nSpikes'),1)>=100,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

function [group,groupid]=interneuronfilter(group,groupid,varnames)
% Only keep putative interneurons 
%
% 1.peak to valley duration above .2ms
% 2.thetaindex about 2
idx=sum(group(:,contains(varnames,'spikewidth'),1)<.2,3)>0 &...
    sum(group(:,contains(varnames,'thetaindex'),1)>=2,3)>0 &...
    sum(group(:,contains(varnames,'nSpikes'),1)>200,3)>0;

groupid=groupid(idx,:);
group=group(idx,:,:);
end

function visualizecells(groupid,group)
cd('F:\ClarkP30_Recordings\ProcessedData')
for i=1:length(groupid)
    data=load(groupid{i,1});
%     close all
    p=postprocessFigures(data,{groupid{i,2},str2double(groupid(i,3))});
    set(gcf, 'Position', get(0, 'Screensize'));
%     pause(3)


    print(gcf,'-dpng', '-r300',...
        ['d:\Users\BClarkLab\Desktop\Laura Temp\PlaceCells\',group,...
        filesep,groupid{i,1},groupid{i,2},groupid{i,3},'.png'])
    close all
    


%     saveas(gcf,['D:\Projects\PAE_PlaceCell',filesep,'test.emf'])
%         p=postprocessFigures2(data,{groupid{i,2},str2double(groupid(i,3))});

% print('-dmeta',['D:\Projects\PAE_PlaceCell',filesep,'test.emf'])
end
end


function autoCorrMat=popVecAutocorr(groupid)
cd('D:\ClarkP30_Recordings\ProcessedData')
autoCorrMat=[];
for i=1:length(groupid)
    data=load(groupid{i,1},'thetaautocorr','spikesID');
    %     close all
    cells=find(contains(data.spikesID.TetrodeNum,groupid{i,2}) & ismember(data.spikesID.CellNum,str2double(groupid(i,3))))';
    
    for ii=cells
        autoCorrMatTEMP=data.thetaautocorr{cells,1}';
        autoCorrMat=[autoCorrMat; autoCorrMatTEMP];  
    end
end

end

function tuning=get_tuning(groupid)
cd('D:\ClarkP30_Recordings\ProcessedData')

sessions={'S1','S2','S3','S4'};
for i=1:length(groupid)
    data=load(groupid{i,1},'events','spikesID','Spikes','frames','samplerate');
    cells=find(contains(data.spikesID.TetrodeNum,groupid{i,2}) & ismember(data.spikesID.CellNum,str2double(groupid(i,3))))';
    
    if size(data.events,2)<4
        continue
    end
    % GET TUNNING CURVE
    for ii=cells
        for iii=1:4
            [data_video_spk,~]=createframes_w_spikebinary(data,iii,ii);
            [~,~,~,~,~,temptuning]=tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),...
                data_video_spk(data_video_spk(:,6)==1,4),data.samplerate);
            if exist('tuning','var') && isfield(tuning,(sessions{iii})) 
                r=size(tuning.(sessions{iii}),1)+1;
            else
                r=1;
            end
            tuning.(sessions{iii})(r,:)=rescale(temptuning,0,1);
        end 
    end
   
 
%     polarplot(deg2rad(0:6:360),tuning,'k')
%     ax=gca;
%     set(ax,'RGrid','off','ThetaGrid','off','ThetaTick',[0 90],...
%         'ThetaTickLabels',[0 90],'RTick',max(tuning),'RTickLabel',max(tuning),'LineWidth',3);
%     axis tight
end

end

% FOR DIAG POP VECTOR
function [mat]=arrangenorm(mat)
fields=fieldnames(mat);
[~,I]=max(mat.(fields{1}),[],2);
[~,I2]=sort(I);
for i=1:length(fields)
    mat.(fields{i})=mat.(fields{i})(I2,:);
end
end
       