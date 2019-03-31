% analyze_measures

% load('/Users/RyanHarvey/Downloads/data (3).mat')
% close all
control={'RH13','RH14','LS21','LS23','LE2821','LE2823'};
PAE={'RH11','RH16','LS17','LS19','LE2813'};

% NUMBER OF CLUSTERS
[nclustersC,nsessionsC,ncellvecC,idsessionC]=countclusters(control,data);
[nclustersP,nsessionsP,ncellvecP,idsessionP]=countclusters(PAE,data);
disp(['Session Count: Control: ',num2str(nsessionsC),', PAE: ',num2str(nsessionsP)])
disp(['Cluster Count: Control: ',num2str(nclustersC),', PAE: ',num2str(nclustersP)])

% EXTRACT MEASURES FOR LINEAR TRACK
[group1,movementc,idsc]=extractLINEAR(control,data);
[group2,movementp,idsp]=extractLINEAR(PAE,data);

% EXTRACT MEASURES FOR CYLINDER
[group1cylinder,idsccylinder]=extractCYLINDER(control,data);
[group2cylinder,idspcylinder]=extractCYLINDER(PAE,data);

% FILTER OUT BAD QUALITY & LOW SPIKING CELLS
[group1,idsc]=qualityfilter(group1,idsc);
[group2,idsp]=qualityfilter(group2,idsp);
[group1cylinder,idsccylinder]=qualityfilter(group1cylinder,idsccylinder);
[group2cylinder,idspcylinder]=qualityfilter(group2cylinder,idspcylinder);


% FILTER OUT INTERNEURONS
% [group1,idsc]=interneuronfilter(group1,idsc);
% [group2,idsp]=interneuronfilter(group2,idsp);
% [group1cylinder,idsccylinder]=interneuronfilter(group1cylinder,idsccylinder);
% [group2cylinder,idspcylinder]=interneuronfilter(group2cylinder,idspcylinder);

%% OPTION TO PLOT SPIKES ON PATH, RATE MAPS, AND PHASE PRECESSION 
% plotphase(data,group1,idsc)
% plotphase(data,group1,idsc)

%%

% FILTER OUT NON-PLACE CELLS
[group1,idsc]=placefilter(group1,idsc);
[group2,idsp]=placefilter(group2,idsp);
[group1cylinder,idsccylinder]=placefilter(group1cylinder,idsccylinder);
[group2cylinder,idspcylinder]=placefilter(group2cylinder,idspcylinder);

% FILTER OUT NAN COLUMNS FOR LINEAR TRACK
varnames=data.RH11.S20160511134416.varnames(:)';
varnames(:,sum(isnan(group1))==size(group1,1))=[];
group1(:,sum(isnan(group1))==size(group1,1))=[];
group2(:,sum(isnan(group2))==size(group2,1))=[];

group1(isinf(group1))=0;
group2(isinf(group2))=0;

% DIRECTIONALITY HISTOGRAM
curfig=figure;
curfig.Color=[1 1 1];
[h1,h2]=CoolHistogram(group1(:,18),group2(:,18),30,varnames{18});


AllStats=CDFplots(group1,group2,{'Sacc','PAE'},varnames,1)


savefigs=0;
if savefigs==1
    for i=1:length(varnames)
        AllStats=CDFplots(group1(:,i),group2(:,i),{'Sacc','PAE'},varnames(i),2);
        curfig=figure(1);
        curfig.Color=[1 1 1];
        %     curfig.OuterPosition=[1 1 960 1058]
        curfig.OuterPosition=[817 404 554 555];
        ax=gca;
        %     set(ax,'Position',[0.13 0.11 0.775 0.815])
        print(curfig,'-bestfit', '-dpdf', '-r300',...
            ['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_retreat/linearTrackMeasures',...
            filesep,[varnames{i},'_fig.pdf']])
        close all
    end
end
AllStats=ScatterBox(movementc(:,1),movementp(:,1),{'Sacc','PAE'},{'Average Velocity (cm/sec)'},2)
% print(figure(1),'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub',filesep,'Average velocity_fig.pdf'])


% FILTER OUT NAN COLUMNS FOR LINEAR CYLINDER
varnames=data.RH11.S20160511134416.varnames(:)';
varnames(:,sum(isnan(group1cylinder(:,:,1)))==size(group1cylinder,1))=[];
group1cylinder(:,sum(isnan(group1cylinder(:,:,1)))==size(group1cylinder,1),:)=[];
group2cylinder(:,sum(isnan(group2cylinder(:,:,1)))==size(group2cylinder,1),:)=[];

% absolute value of phase corr
% group1cylinder(:,13,:)=abs(group1cylinder(:,12,:));
% group2cylinder(:,13,:)=abs(group2cylinder(:,12,:));
group1cylinder(isinf(group1cylinder))=0;
group2cylinder(isinf(group2cylinder))=0;
group1cylinder(isnan(group1cylinder))=0;
group2cylinder(isnan(group2cylinder))=0;

AllStats=CDFplots(group1cylinder(:,:,1),group2cylinder(:,:,1),{'Sacc','PAE'},varnames,1)


%%
% plotphasecirc(data,group1cylinder,idsccylinder)
% plotphasecirc(data,group2cylinder,idspcylinder)

%%
savefigs=0;
if savefigs==1
    for i=1:length(varnames)
        AllStats=CDFplots(group1cylinder(:,i),group2cylinder(:,i),{'Sacc','PAE'},varnames(i),2);
        curfig=figure(1);
        curfig.Color=[1 1 1];
        %     curfig.OuterPosition=[1 1 960 1058]
        curfig.OuterPosition=[817 404 554 555];
        ax=gca;
        %     set(ax,'Position',[0.13 0.11 0.775 0.815])
        print(curfig,'-bestfit', '-dpdf', '-r600',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/Learning&Memory/ArenaMeasures',filesep,[varnames{i},'_fig.pdf']])
        close all
    end
end

% \\\\\\\\\\\\\\\\\\\\ LOCAL FUNCTIONS BELOW \\\\\\\\\\\\\\\\\\\\

function [cells,movement,ids]=extractLINEAR(rats,data)
cells=[];
movement=[];
ids=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        ncells=fieldnames(data.(rats{r}).(sessions{s}).('thetaautocorr'));
        ids=[ids;repmat({(rats{r}),(sessions{s})},length(ncells),1),ncells,repmat({'session1'},length(ncells),1)];
        ids=[ids;repmat({(rats{r}),(sessions{s})},length(ncells),1),ncells,repmat({'session2'},length(ncells),1)];
        movement=[movement;data.(rats{r}).(sessions{s}).BasicLoco.AverageAnglePerSec,data.(rats{r}).(sessions{s}).BasicLoco.OverallDistTraveled,data.(rats{r}).(sessions{s}).BasicLoco.MeanVelocity];
        cells=[cells;[data.(rats{r}).(sessions{s}).measures(:,:,1);data.(rats{r}).(sessions{s}).measures(:,:,2)]];
    end
end
end

function [cells,ids]=extractCYLINDER(rats,data)
cells=[];
ids=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        [~,~,d]=size(data.(rats{r}).(sessions{s}).measures);
        if d>3
            ncells=fieldnames(data.(rats{r}).(sessions{s}).('thetaautocorr'));
            ids=[ids;repmat({(rats{r}),(sessions{s})},length(ncells),1),ncells,repmat({'session3'},length(ncells),1)];
            ids=[ids;repmat({(rats{r}),(sessions{s})},length(ncells),1),ncells,repmat({'session4'},length(ncells),1)];
            cells=[cells;data.(rats{r}).(sessions{s}).measures(:,:,[3:4])];
        end
    end
end
end

function [groupout,ids]=placefilter(groupin,ids)
% Only keep cells with the following
%
% 1. > .25 info content
% 2. > .3 average rate
% 3. > 2hz peak rate

groupout=groupin(groupin(:,1,1)>.25 & groupin(:,5,1)>.3 & groupin(:,4,1)>2,:,:);
ids=ids(groupin(:,1,1)>.25 & groupin(:,5,1)>.3 & groupin(:,4,1)>2,:,:);
end

function [groupout,ids]=qualityfilter(groupin,ids)
% Only keep cells with the following characteristics
%
% 1. above 50 spikes
% 2. below 0.5% of spikes in the first 2ms on the autocorrelation
% 3. below 50% cut off by threshold on plots of peak spike amplitude on one tetrode wire versus another
% 4. filter linear track by at least 5 laps
% 5. > .02 temporal stability (checks for single burst artifacts... 0.02 is a light cut off)


% groupout=groupin(groupin(:,8,1)>50 & groupin(:,43,1)<.05 & groupin(:,44,1)<.5 & groupin(:,59,1)>5 | isnan(groupin(:,59,1)) & groupin(:,58,1)>.02,:,:);
% ids=ids(groupin(:,8,1)>50 & groupin(:,43,1)<.05 & groupin(:,44,1)<.5 & groupin(:,59,1)>5 | isnan(groupin(:,59,1)) & groupin(:,58,1)>.02,:,:);

groupout=groupin(groupin(:,8,1)>50 & [groupin(:,58,1)>5 | isnan(groupin(:,58,1))] & groupin(:,28,1)>0,:,:);
ids=ids(groupin(:,8,1)>50 & [groupin(:,58,1)>5 | isnan(groupin(:,58,1))] & groupin(:,28,1)>0,:,:);

end

function [groupout,ids]=interneuronfilter(groupin,ids)
% Only keep pyramidal cells
%
% 1. average rate below 10hz
% 2. peak to valley duration above .2ms
groupout=groupin(groupin(:,5,1)<10 & groupin(:,52,1)>0.2,:,:);
ids=ids(groupin(:,5,1)<10 & groupin(:,52,1)>0.2,:,:);

end

function [ncellz,nsessions,ncells,ids]=countclusters(rats,data)
ncells=[];
ids=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        ids=[ids;{(rats{r}),(sessions{s})}];
        ncells=[ncells;length((data.(rats{r}).(sessions{s}).Spikes))];
    end
end
nsessions=length(ncells);
ncellz=sum(ncells);
end
% function [cells,movement,ids]=findprinciplecells(rats,data)
% cells=[];
% movement=[];
% ids=[];
% for r=1:length(rats)
%     sessions=fieldnames(data.(rats{r}));
%     for s=1:length(sessions)
%
%
%
%         for ch=1:length(data.(rats{r}).(sessions{s}).avgwave)
%             tetrode=data.(rats{r}).(sessions{s}).avgwave{ch};
%             max(tetrode')
%         end
%
%         ids=[ids;repmat({(rats{r}),(sessions{s})},length(ncells),1),ncells,repmat({'session1'},length(ncells),1)];
%         ids=[ids;repmat({(rats{r}),(sessions{s})},length(ncells),1),ncells,repmat({'session2'},length(ncells),1)];
%         movement=[movement;data.(rats{r}).(sessions{s}).BasicLoco.AverageAnglePerSec,data.(rats{r}).(sessions{s}).BasicLoco.OverallDistTraveled,data.(rats{r}).(sessions{s}).BasicLoco.MeanVelocity];
%         cells=[cells;[data.(rats{r}).(sessions{s}).measures(:,:,1);data.(rats{r}).(sessions{s}).measures(:,:,2)]];
%     end
% end
% end
function plotphase(data,group,ids)

close all
I=find(group(:,1)>.2 & group(:,4)>1.5 & group(:,7)<30 & group(:,35)>50);


potentialcells=ids(I,:);

for i=1:length(potentialcells)
    
    % EXTRACT FRAME AND SPIKE DATA TO PLOT TRACK
    frames=data.(potentialcells{i,1}).(potentialcells{i,2}).frames;
    % get spikes
    cellnum=regexp(potentialcells(i,3),'\d*','Match');
    spkts=data.(potentialcells{i,1}).(potentialcells{i,2}).Spikes{str2double(cellnum{1})};
    % get events
    events=data.(potentialcells{i,1}).(potentialcells{i,2}).events;
    
    % find track length
    date=char(potentialcells{i,2});
    date=[date(2:5),'-',date(6:7),'-',date(8:9),'_',date(10:11),'-',date(12:13),'-',date(14:15)];
    track_length=TrackLength(['/Volumes/Ryan_4TB/Place_Cell_Data/RawPAE_PlaceCell',...
        filesep,potentialcells{i,1},filesep,date]);
    
    [right,left]=rebuildFramematrix(frames,spkts,events,track_length);
    clear date events spkts frames
    
    
    if length(data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),1})>2
        xl=data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),1}(:,1);
        pl=data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),1}(:,2);
    else
        xl=NaN;
        pl=NaN;
    end
    
    if length(data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),2})>2
        xr=data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),2}(:,1);
        pr=data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),2}(:,2);
    else
        xr=NaN;
        pr=NaN;
    end
    
    ratemap1=data.(potentialcells{i,1}).(potentialcells{i,2}).ratemap.(char(strcat('Cell',cellnum{1}))).session1;
    ratemap2=data.(potentialcells{i,1}).(potentialcells{i,2}).ratemap.(char(strcat('Cell',cellnum{1}))).session2;
    PR=max([ratemap1,ratemap2]);
    
    fig=figure('Name',strcat(potentialcells{i,:}));fig.Color=[1 1 1];fig.OuterPosition=[1 6 1117 628];
    subplot(4,2,1)
    plot(right.dataspks(:,2),right.dataspks(:,1),'.k');hold on
    scatter(right.dataspks(right.dataspks(:,6)==1 & right.dataspks(:,5)>5,2),right.dataspks(right.dataspks(:,6)==1 & right.dataspks(:,5)>5,1),'filled','r')
    ylabel('Laps')
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTick',...
        linspace(min(right.dataspks(:,2)),max(right.dataspks(:,2)),3),...
        'XTickLabel',[],'box','off','YTickLabel',[],'TickLength',[0 0])
    xlim([min(right.dataspks(:,2)),max(right.dataspks(:,2))])
    ylim([min(right.dataspks(:,1)),max(right.dataspks(:,1))])
    
    subplot(4,2,2)
    plot(left.dataspks(:,2),left.dataspks(:,1),'.k');hold on
    scatter(left.dataspks(left.dataspks(:,6)==1 & left.dataspks(:,5)>5,2),left.dataspks(left.dataspks(:,6)==1 & left.dataspks(:,5)>5,1),'filled','r')
    ylabel('Laps')
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTick',...
        linspace(min(left.dataspks(:,2)),max(left.dataspks(:,2)),3),...
        'XTickLabel',[],'box','off','YTickLabel',[],'TickLength',[0 0])
    xlim([min(left.dataspks(:,2)),max(left.dataspks(:,2))])
    ylim([min(left.dataspks(:,1)),max(left.dataspks(:,1))])
    
    subplot(4,2,3)
    rainbowplot(ratemap1);
    ylim([0 PR])
    ylabel('FR (hz)')
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTickLabel',[],'box','off')
    
    subplot(4,2,4)
    rainbowplot(ratemap2);
    ylim([0 PR])
    ylabel('FR (hz)')
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTickLabel',[],'box','off')
    
    subplot(4,2,5)
    scatter(xl,pl,'filled','k');hold on
    scatter(xl,pl+360,'filled','k')
    xlim([0 1])
    ylim([0 720])
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTickLabel',[],'box','off',...
        'YTick',[0 360 720],'YTickLabel',[0 360 720],'TickLength',[0 0])
    
    subplot(4,2,6)
    scatter(xr,pr,'filled','k');hold on
    scatter(xr,pr+360,'filled','k')
    xlim([0 1])
    ylim([0 720])
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTickLabel',[],'box','off',...
        'YTick',[0 360 720],'YTickLabel',[0 360 720],'TickLength',[0 0])
    
    
    filtWidth = [360 7]; filtSigma = 40;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    
    subplot(4,2,7)
    [N,~,~]=histcounts2([pl;pl+360],[xl;xl],linspace(0,720,720+1),linspace(0,1,100+1));
    N = nanconv([N;N;N],imageFilter, 'nanout');
    N=N(721:1440,:);
    pcolor(N);
    shading interp
    axis xy
    colormap jet
    xlabel('Distance on track (cm)')
    ylabel('Theta Phase (deg)')
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTick',...
        [1 50.5 100],'XTickLabel',[0 60 120],'box','off',...
        'YTick',[1 360.5 720],'YTickLabel',[0 360 720],'TickLength',[0 0])
    
    subplot(4,2,8)
    [N,~,~]=histcounts2([pr;pr+360],[xr;xr],linspace(0,720,720+1),linspace(0,1,100+1));
    N = nanconv([N;N;N],imageFilter, 'nanout');
    N=N(721:1440,:);
    pcolor(N);
    shading interp
    axis xy
    colormap jet
    xlabel('Distance on track (cm)')
    ylabel('Theta Phase (deg)')
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTick',...
        [1 50.5 100],'XTickLabel',[0 60 120],'box','off',...
        'YTick',[1 360.5 720],'YTickLabel',[0 360 720],'TickLength',[0 0])
    
    print(fig,'-dpng', '-r300',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_retreat/Figures/Control',filesep,potentialcells{i,:},'.png'])

end
end

function [right,left]=rebuildFramematrix(frames,SpikeFile,events,track_length)
% restrict to linear track session
if sum(events==[1;1])~=2 % if more than one session
    frames=frames(frames(:,1)>events(1,1) & frames(:,1)<events(2,1),:);
end

% calc head angle
ExtractedAngle=XYangle(frames(:,2),frames(:,3));

% calc vel
vel_cmPerSec=abs(diff(frames(:,2)))*track_length/(range(frames(:,2)))*30;

frames(:,[4:5])=[[ExtractedAngle;ExtractedAngle(end)],[vel_cmPerSec;vel_cmPerSec(end)]];

if sum(events==[1;1])~=2 % if more than one session
    SpikeFile=SpikeFile(SpikeFile(:,1)>events(1,1) & SpikeFile(:,1)<events(2,1),:);
end

% in=contiguousframes(frames(:,5)<2,6);
data_video_nospk=frames;

% INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
% TS=interp1(data_video_nospk(:,1),data_video_nospk(:,1),SpikeFile,'linear');
X=interp1(data_video_nospk(:,1),data_video_nospk(:,2),SpikeFile,'linear');
Y=interp1(data_video_nospk(:,1),data_video_nospk(:,3),SpikeFile,'linear');
A=interp1(data_video_nospk(:,1),data_video_nospk(:,4),SpikeFile,'linear');
VEL=interp1(data_video_nospk(:,1),data_video_nospk(:,5),SpikeFile,'linear');
% VELidx=interp1(data_video_nospk(:,1),data_video_nospk(:,6),SpikeFile,'linear');

% CONCAT AND SORT
data_video_spk=sortrows([[SpikeFile X Y A VEL ones(size(SpikeFile,1),1)];[data_video_nospk,zeros(length(data_video_nospk),1)]],1);


% SPLIT RIGHT AND LEFT RUNS
[right,left,~,~,~]=RightVsLeft(data_video_nospk,data_video_spk,track_length,30);


% [lapsR,mapsR,corsR]=splitlapagain(right.dataspks,find(diff(right.dataspks(:,1))>1e7),track_length,right.SmoothRateMap);
% [lapsL,mapsL,corsL]=splitlapagain(left.dataspks,find(diff(left.dataspks(:,1))>1e7),track_length,left.SmoothRateMap);
end

function plotphasecirc(data,group,ids)

close all
% I=find(group(:,1,1)>.2 & group(:,4,1)>1.5 & group(:,7,1)<30 & group(:,35,1)>50);
% potentialcells=ids(I,:);

potentialcells=ids;

for i=1:length(potentialcells)

    cellnum=regexp(potentialcells(i,3),'\d*','Match');

    if length(data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),3})>2
        x=data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),3}(:,1);
        p=data.(potentialcells{i,1}).(potentialcells{i,2}).ThPrecess{str2double(cellnum{1}),3}(:,2);
    else
        x=NaN;
        p=NaN;
    end
    
    fig=figure('Name',strcat(potentialcells{i,:}));fig.Color=[1 1 1];fig.OuterPosition=[1 6 1117 628];

    subplot(2,1,1)
    scatter(x,p,'filled','k');hold on
    scatter(x,p+360,'filled','k')
    xlim([0 1])
    ylim([0 720])
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTickLabel',[],'box','off',...
        'YTick',[0 360 720],'YTickLabel',[0 360 720],'TickLength',[0 0])
   
    filtWidth = [360 7]; filtSigma = 40;
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    
    subplot(2,1,2)
    [N,~,~]=histcounts2([p;p+360],[x;x],linspace(0,720,720+1),linspace(0,1,100+1));
    N = nanconv([N;N;N],imageFilter, 'nanout');
    N=N(721:1440,:);
    pcolor(N);
    shading interp
    axis xy
    colormap jet
    xlabel('Distance')
    ylabel('Theta Phase (deg)')
    set(gca,'FontWeight','bold','FontSize',12,'LineWidth',3,'XTick',...
        [1 50.5 100],'XTickLabel',[0 .5 1],'box','off',...
        'YTick',[1 360.5 720],'YTickLabel',[0 360 720],'TickLength',[0 0])

%     print(fig,'-dpng', '-r300',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_retreat/Figures/Control',filesep,potentialcells{i,:},'.png'])

end
end


