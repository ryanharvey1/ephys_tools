% analyze_measures

load('F:\P30_Recordings/data.mat')
close all
control={'LB01','LB03','LB05'};
TG={'LB02','LB04','LB06'};

[shuffDist]=extractSHUFF([control,TG],data);
cutoff=prctile(shuffDist,95);

% EXTRACT MEASURES FOR CYLINDER
[WTrats,idsccylinder]=extractCYLINDER(control,data);
[TGrats,idspcylinder]=extractCYLINDER(TG,data);

% FILTER OUT BAD QUALITY & LOW SPIKING CELLS
[WTrats,idsccylinder]=qualityfilter(WTrats,idsccylinder);
[TGrats,idspcylinder]=qualityfilter(TGrats,idspcylinder);

% 
% % FILTER OUT INTERNEURONS
% [WTrats,idsccylinder]=interneuronfilter(WTrats,idsccylinder);
% [TGrats,idspcylinder]=interneuronfilter(TGrats,idspcylinder);

% % CHECK OUT PLACE MODULATED CELLS
% [WTrats,idsccylinder]=Placefilter(WTrats,idsccylinder);
% [TGrats,idspcylinder]=Placefilter(TGrats,idspcylinder);


% FILTER OUT NON-DIRECTIONAL CELLS
[WTrats,idsccylinder]=HDfilter(WTrats,idsccylinder,cutoff);
[TGrats,idspcylinder]=HDfilter(TGrats,idspcylinder,cutoff);
% 
% plotHD(idsccylinder,WTrats,data)
plotHD(idspcylinder,TGrats,data)

% FILTER OUT NAN COLUMNS FOR LINEAR CYLINDER
varnames=data.LB03.S20170309182117.varnames(:)';
varnames(:,sum(isnan(WTrats(:,:,1)))==size(WTrats,1))=[];
WTrats(:,sum(isnan(WTrats(:,:,1)))==size(WTrats,1),:)=[];
TGrats(:,sum(isnan(TGrats(:,:,1)))==size(TGrats,1),:)=[];

% absolute value of phase corr
% WTrats(:,13,:)=abs(WTrats(:,12,:));
% TGrats(:,13,:)=abs(TGrats(:,12,:));
WTrats(isinf(WTrats))=0;
TGrats(isinf(TGrats))=0;
WTrats(isnan(WTrats))=0;
TGrats(isnan(TGrats))=0;

%PLOT STATS
for i=1:4
figure
AllStats=CDFplots(WTrats(:,:,1),TGrats(:,:,1),{'WT','TG'},varnames,2)
end
% \\\\\\\\\\\\\\\\\\\\ LOCAL FUNCTIONS BELOW \\\\\\\\\\\\\\\\\\\\


function [cells,ids]=extractCYLINDER(rats,data)
cells=[];
ids=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        if ~isfield(data.(rats{r}).(sessions{s}),'measures')
            continue
        end
        [row,col,d]=size(data.(rats{r}).(sessions{s}).measures);
  
        tempdat=NaN(row,col,6);
        tempdat(:,:,1:d)= data.(rats{r}).(sessions{s}).measures(:,:,:);
        cells=[cells;tempdat];
        ids=[ids;[repmat({(rats{r}),(sessions{s})},row,1),num2cell([1:row]')]];
    end
end
end



function [shuffDist]=extractSHUFF(rats,data)
shuffDist=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        if ~isfield(data.(rats{r}).(sessions{s}),'rlength_shuff')
            continue
        end
        for cells=1:size(data.(rats{r}).(sessions{s}).rlength_shuff,1)
            if data.(rats{r}).(sessions{s}).measures(cells,7,1)>50 &&...
                    data.(rats{r}).(sessions{s}).measures(cells,11,1)<.05 &&...
                    data.(rats{r}).(sessions{s}).measures(cells,12,1)<.5
                shuffDist=[shuffDist;data.(rats{r}).(sessions{s}).rlength_shuff{cells,1}];
            end
        end
        
    end
end
end

function [groupout,ids]=HDfilter(groupin,ids,cutoff)
% Only keep cells with the following
%
% 1. pass shuffled distribution 95th percentile (rLength)

% groupout=groupin(sum(groupin(:,2,:)==1,3)>0,:,:);
groupout=groupin(sum(groupin(:,2,2)>=cutoff,3)>0,:,:);
ids=ids(sum(groupin(:,2,2)>=cutoff,3)>0,:,:);
end

function [groupout,ids]=Placefilter(groupin,ids)
% Only keep cells with the following
%
% 1. pass shuffled distribution 95th percentile (rLength)

% groupout=groupin(sum(groupin(:,2,:)==1,3)>0,:,:);
groupout=groupin(sum(groupin(:,32,:)>=.5,3)>0,:,:);
ids=ids(sum(groupin(:,32,:)>=.5,3)>0,:,:);
end


function plotHD(ids,group1cylinder,data)
for i=1:length(ids)
    % EXTRACT FRAME AND SPIKE DATA TO PLOT TRACK
    frames=data.(ids{i,1}).(ids{i,2}).frames;
    % get spikes
    spkts=data.(ids{i,1}).(ids{i,2}).Spikes{ids{i,3}};
    % get events
    events=data.(ids{i,1}).(ids{i,2}).events;
    data_video_spk=rebuildFramematrix(frames,spkts,events);
    if data_video_spk==0
       
        disp(ids{i,:})
         test=1
    end
    da=pi/30;
    angBins=da/2:da:2*pi-da/2;
    hdtun=data.(ids{i,1}).(ids{i,2}).hdTuning{ids{i,3}};
%     rlength=data.(ids{i,1}).(ids{i,2}).measures{ids{i,2}};
%     preferred_Direction=data.(ids{i,1}).(ids{i,2}).measures{ids{i,8}};
    % Firing rate x HD polar plot for the nonsmoothed data above
        figure(i); polarplot = polar(angBins',hdtun','b');
        set(polarplot, 'linewidth',3,'color','k'); axis off
%         title(['Polor Plot, MeanVecLength: ',num2str(rlength),' Pref_Dir: ',num2str(preferred_Direction)]);
        set(0,'Showhiddenhandles','on')
        
        % ---------CODE FOR PUBLICATION FIGURE--------
        extrastuff = setdiff(get(gca,'children'),polarplot);
        delete(extrastuff)
        horizontal=line([-max(hdtun) max(hdtun)],[0 0]);
        vertical=line([0 0],[-max(hdtun) max(hdtun)]);
        set(horizontal,'linewidth',2,'color','k');
        set(vertical,'linewidth',2,'color','k');
        %---------------------------------------------
        set(figure(i),'Position',[686 325 977 619]);
        fig1 = figure(i);
        
        clear angBins
%     
%     tuningfig=figure; tuningfig.Color=[1 1 1];
%     subplot(1,2,1)
%     p=plot(hdtun,'k');hold on
%     xlabel('Head Angle')
%     ylabel('Firing Rate (hz)')
%     title([ids{i,1,1}])
%     set(p,'LineWidth',3)
%     set(gca,'box','off','LineWidth',2,'XTick',linspace(0,60,7),'XTickLabel',linspace(0,60,7)*6,'FontSize',20,'FontWeight','bold')
%     annotation('textbox',...
%         [0.15 0.65 0.3 0.15],...
%         'String',{['r=' num2str(group1cylinder(i,2,1))],...
%         ['hlfWdth=' num2str(group1cylinder(i,3,1))], ['DIC=' num2str(group1cylinder(i,4,1))]...
%         ['nSpike=' num2str(group1cylinder(i,7,1))], ['peakRate=' num2str(group1cylinder(i,5,1))]},...
%         'FontSize',14,...
%         'FontName','helvetica','EdgeColor','none');
%     subplot(1,2,2)
%     plot(data_video_spk(:,2),data_video_spk(:,3),'.k');hold on
%     scatter(data_video_spk(data_video_spk(:,4)==1,2),data_video_spk(data_video_spk(:,4)==1,3),'Filled','r');
%     axis image
%     axis off
    
end
end

function [groupout,ids]=qualityfilter(groupin,ids)
% Only keep cells with the following characteristics
%
% 1. above 50 spikes
% 2. below 0.5% of spikes in the first 2ms on the autocorrelation
% 3. below 50% cut off by threshold on plots of peak spike amplitude on one tetrode wire versus another
%
% before=length(groupin)/2;
groupout=groupin(sum(groupin(:,7,:)>50,3) & sum(groupin(:,11,:)<.05,3) & sum(groupin(:,12,:)<.75,3),:,:);
ids=ids(sum(groupin(:,7,:)>50,3) & sum(groupin(:,11,:)<.05,3) & sum(groupin(:,12,:)<.75,3),:,:);
% after=length(groupout)/2;

% disp(['removing ',num2str(before-after),' bad quality cells'])
end

function [groupout,ids]=interneuronfilter(groupin,ids)
% Only keep pyramidal cells
%
% 1. average rate below 10hz
% 2. peak to valley duration above .2ms
groupout=groupin(groupin(:,6,1)<10 & groupin(:,21,1)>0.2,:,:);
ids=ids(groupin(:,6,1)<10 & groupin(:,21,1)>0.2,:,:);
% groupout=groupin(groupin(:,6,1)>10 & groupin(:,21,1)<0.2,:,:);
% ids=ids(groupin(:,6,1)>10 & groupin(:,21,1)<0.2,:,:);
end

function [data_video_spk]=rebuildFramematrix(frames,SpikeFile,events)
% restrict to linear track session
if sum(events==[1;1])~=2 % if more than one session
    frames=frames(frames(:,1)>events(1,1) & frames(:,1)<events(2,1),:);
end

if sum(events==[1;1])~=2 % if more than one session
    SpikeFile=SpikeFile(SpikeFile(:,1)>events(1,1) & SpikeFile(:,1)<events(2,1),:);
end
if (frames/30)<60
    test=1;
    data_video_spk=0;
    return
end

% INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
TS=interp1(frames(:,1),frames(:,1),SpikeFile,'linear');
X=interp1(frames(:,1),frames(:,2),SpikeFile,'linear');
Y=interp1(frames(:,1),frames(:,3),SpikeFile,'linear');

% CONCAT AND SORT
data_video_spk=sortrows([[TS X Y ones(size(TS,1),1)];[frames,zeros(length(frames),1)]],1);
end

