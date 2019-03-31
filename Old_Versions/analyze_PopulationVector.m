% analyze_PopulationVector
% Directional Population vector analysis for newer SNAP sort data

% load('data.mat')
% addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/gif'))
% close all
control={'RH13','RH14','LS21','LS23','LE2821','LE2823'};
PAE={'RH11','RH16','LS17','LS19','LE2813'};


% EXTRACT RATEMAPS FROM DATA FILE, NORM, AND CORR
[mapright_1C,mapleft_1C,RowCorr_1C,ratediff1C,unnormright1C,unnormleft1C,ids1C]=extract(control,data,1);
[mapright_2C,mapleft_2C,RowCorr_2C,ratediff2C,unnormright2C,unnormleft2C,ids2C]=extract(control,data,2);
[mapright_1P,mapleft_1P,RowCorr_1P,ratediff1P,unnormright1P,unnormleft1P,ids1P]=extract(PAE,data,1);
[mapright_2P,mapleft_2P,RowCorr_2P,ratediff2P,unnormright2P,unnormleft2P,ids2P]=extract(PAE,data,2);

% collapse ids so that we are sampling without replacement
clear nids1C nids2C nids1P nids2P
for i=1:length(ids1C)
    nids1C{i,1}=strcat(ids1C{i,:});
end
for i=1:length(ids2C)
    nids2C{i,1}=strcat(ids2C{i,:});
end
for i=1:length(ids1P)
    nids1P{i,1}=strcat(ids1P{i,:});
end
for i=1:length(ids2P)
    nids2P{i,1}=strcat(ids2P{i,:});
end


% % WORKING
% placemaps_no_repeatsC=[mapright_1C(~ismember(nids1C,nids2C),:);mapright_2C];
% oppositemaps_no_repeatsC=[mapleft_1C(~ismember(nids1C,nids2C),:);mapleft_2C];
% 
% placemaps_no_repeatsP=[mapright_1P(~ismember(nids1P,nids2P),:);mapright_2P];
% oppositemaps_no_repeatsP=[mapleft_1P(~ismember(nids1P,nids2P),:);mapleft_2P];

% for i=1:length(placemaps_no_repeatsC)
% [m,I(i,1)]=max(placemaps_no_repeatsC(i,:));
% end
% histogram(I)
% 
% for i=1:length(placemaps_no_repeatsC)
% [m,I(i,1)]=max(oppositemaps_no_repeatsC(i,:));
% end
% histogram(I)
% 
% 
% for i=1:length(placemaps_no_repeatsP)
% [m,I(i,1)]=max(placemaps_no_repeatsP(i,:));
% end
% 
% for i=1:length(placemaps_no_repeatsP)
% [m,I(i,1)]=max(oppositemaps_no_repeatsP(i,:));
% end

%% plot place fields and their opposite runs
[allCR,allCL]=arrange([mapright_1C;mapleft_2C],[mapleft_1C;mapright_2C]);
[allPR,allPL]=arrange([mapright_1P;mapleft_2P],[mapleft_1P;mapright_2P]);
subplot(2,2,1)
imagesc([allCR,allCL]);colormap jet
subplot(2,2,2)
imagesc([allPR,allPL]);colormap jet
subplot(2,2,3)
imagesc(corr(allCR,allCL));colormap jet
subplot(2,2,4)
imagesc(corr(allPR,allPL));colormap jet

figure;imagesc(corr(allCR,allCR));colormap jet
figure;imagesc(corr(allPR,allPR));colormap jet



% subplot(1,2,1)
% imagesc([allCR,alCL]);colormap jet;set(gca,'Visible','off','box','off');
% subplot(1,2,2)
% imagesc([allPR,allPL]);colormap jet;set(gca,'Visible','off','box','off');

%% REARANGE RATEMAPS IN DIAG
[mapright_1C,mapleft_1C]=arrange(mapright_1C,mapleft_1C);
[mapleft_2C,mapright_2C]=arrange(mapleft_2C,mapright_2C);
[mapright_1P,mapleft_1P]=arrange(mapright_1P,mapleft_1P);
[mapleft_2P,mapright_2P]=arrange(mapleft_2P,mapright_2P);

%% SEE IF FIELDS ARE REMAPPING BASED ON THEIR DISTANCE FROM THE END OF THE
% TRACK i.e. if they are using local cues to stabilze their fields
% flipcor1C=flipcorr(mapright_1C,mapleft_1C);
% flipcor2C=flipcorr(mapleft_2C,mapright_2C);
% flipcor1P=flipcorr(mapright_1P,mapleft_1P);
% flipcor2P=flipcorr(mapleft_2P,mapright_2P);
%
% flipmapstat=CDFplots([flipcor1C;flipcor2C],[flipcor1P;flipcor2P],{'Sacc','PAE'},{'Flip Map Corr(r)'},1)
% normalmapcorstats=CDFplots([RowCorr_1C;RowCorr_2C],[RowCorr_1P;RowCorr_2P],{'Sacc','PAE'},{'Directional Corr(r)'},1)
%


%% RATEOVERLAP/FIELDOVERLAP
[fieldoverlapC,DisplacementC,rateoverlapC,distalcodeC,positioncodeC]=OVERLAP([ids1C(~ismember(nids1C,nids2C),:);ids2C],data);
[fieldoverlapP,DisplacementP,rateoverlapP,distalcodeP,positioncodeP]=OVERLAP([ids1P(~ismember(nids1P,nids2P),:);ids2P],data);
%%
% CDFplots(fieldoverlapC,fieldoverlapP,{'Sacc','PAE'},{'Field Overlap'},1)
% CDFplots(DisplacementC,DisplacementP,{'Sacc','PAE'},{'Displacement (cm)'},1)
% CDFplots(rateoverlapC,rateoverlapP,{'Sacc','PAE'},{'Rate Overlap'},2)
label=0;
if label==1
labels={'justfieldC','justrateC','justdistalC','justpositionC','ratebyfieldC','ratebydistalC','unidentified'};
[idxmatC]=labeldirectionalcodes(positioncodeC,distalcodeC,rateoverlapC);
[idxmatP]=labeldirectionalcodes(positioncodeP,distalcodeP,rateoverlapP);
%
for i=1:size(idxmatC,2)
    [h,p, chi2stat,df]=prop_test([sum(idxmatC(:,i)),sum(idxmatP(:,i))],[length(idxmatC),length(idxmatP)],0);
    disp([labels{i},' X2(',num2str(df),')= ',num2str(chi2stat),' p= ',num2str(p)])
end

piefig=figure; piefig.Color=[1 1 1];
ax=subplot(1,2,1);
p1=pie(sum(idxmatC,1));
title(ax,'Control');ax.FontSize=20;
hold on
ax=subplot(1,2,2);
p2=pie(sum(idxmatP,1));
title(ax,'PAE');ax.FontSize=20;
% print(figure(3),'-bestfit', '-dpdf', '-r600',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/Learning&Memory/Directionality_Figures',filesep,'PIE_fig.pdf'])


histfig=figure;histfig.Color=[1 1 1];
h1=histogram(fieldoverlapC,20,'Normalization','probability');hold on
h2=histogram(fieldoverlapP,20,'Normalization','probability');
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel('Field Overlap')
ylabel('Probability')
end
%%
% controlids=[ids1C(~ismember(nids1C,nids2C),:);ids2C];
% paeids=[ids1P(~ismember(nids1P,nids2P),:);ids2P];
% 
% controlids(idxmatC(:,3)==1,:)
% paeids(idxmatP(:,3)==1,:)
% for i=1:length(controlids)
%     cellnum=regexp(controlids(i,3),'\d*','Match');        
%     ratemap1=data.(controlids{i,1}).(controlids{i,2}).ratemap.(char(strcat('Cell',cellnum{1}))).session1
%     ratemap2=data.(controlids{i,1}).(controlids{i,2}).ratemap.(char(strcat('Cell',cellnum{1}))).session2;
% end



%% FIND BETWEEN LAP STABILITY (the thought is that PAE cells come online later in the
stability=1;
if stability==1
%% session compared with Control cells)

% betweenlapstability([{'LS19','S20170522113749','Cell13'}],data)

[stabilityC,time2stableC,firststablelapC,meanstabilityC,lap_stabilityC]=betweenlapstability([ids1C(~ismember(nids1C,nids2C),:);ids2C],data);
[stabilityP,time2stableP,firststablelapP,meanstabilityP,lap_stabilityP]=betweenlapstability([ids1P(~ismember(nids1P,nids2P),:);ids2P],data);
%%
lapstabilitystat=CDFplots([stabilityC(:,1);stabilityC(:,2)],[stabilityP(:,1);stabilityP(:,2)],{'Sacc','PAE'},{'lap-stability (r)'},2)
time2stablestat=CDFplots([time2stableC(:,1);time2stableC(:,2)],[time2stableP(:,1);time2stableP(:,2)],{'Sacc','PAE'},{'Time to Stability (s)'},2)
firststablelapstat=CDFplots([firststablelapC(:,1);firststablelapC(:,2)],[firststablelapP(:,1);firststablelapP(:,2)],{'Sacc','PAE'},{'lap2stablestat (laps)'},2)

firststablelapstat=CDFplots(meanstabilityC,meanstabilityP,{'Sacc','PAE'},{'Mean Stability'},2)
firststablelapstat=CDFplots(lap_stabilityC,lap_stabilityP,{'Sacc','PAE'},{'Lap Stability'},2)


% STABILITY HISTOGRAM
histfig=figure;histfig.Color=[1 1 1];
h1=histogram([stabilityC(:,1);stabilityC(:,2)],20,'Normalization','probability');hold on
h2=histogram([stabilityP(:,1);stabilityP(:,2)],20,'Normalization','probability');
p1=plot([-.2,-.2],[0,.15],[.2,.2],[0,.15]);hold on
set(p1,'Color','k','LineWidth',3,'LineStyle','--')
set(h1,'FaceColor',[.1 .1 .1],'EdgeColor','k')
set(h2,'FaceColor','r','EdgeColor','k')
set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3)
xlabel('Correlation (r)')
ylabel('Probability')
% print(histfig,'-bestfit', '-dpdf', '-r600',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/Learning&Memory/',filesep,'corrhist_fig.pdf'])

% PORPORTION OF STABILITY
CS=[stabilityC(:,1);stabilityC(:,2)];
PS=[stabilityP(:,1);stabilityP(:,2)];
CS(isnan(CS))=[];
PS(isnan(PS))=[];

n1 = sum(CS<-.2); N1 = length(CS);
n2 = sum(CS>-.2 & CS<.2); N2 = length(CS);
n3 = sum(CS>.2); N3 = length(CS);
n5 = sum(PS<-.2); N5 = length(PS);
n6 = sum(PS>-.2 & PS<.2); N6 = length(PS);
n7 = sum(PS>.2); N7 = length(PS);
x1 = [repmat('a',N1,1); repmat('b',N2,1); repmat('c',N3,1); repmat('e',N5,1); repmat('f',N6,1); repmat('g',N7,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1); repmat(1,n3,1); repmat(2,N3-n3,1);...
    repmat(1,n5,1); repmat(2,N5-n5,1); repmat(1,n6,1); repmat(2,N6-n6,1); repmat(1,n7,1); repmat(2,N7-n7,1)];
[tbl,chi2stat,pval,~] = crosstab(x1,x2)
% below -.2
n1 = sum(CS<-.2); N1 = length(CS);
n2 = sum(PS<-.2); N1 = length(PS);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
% >-.2 & <.2
n1 = sum(CS>-.2 & CS<.2); N2 = length(CS);
n2 = sum(PS>-.2 & PS<.2); N2 = length(PS);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
% >.2
n1 = sum(CS>.2); N1 =  length(CS);
n2 = sum(PS>.2); N2 =  length(PS);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)

clear n1 n2 n3 n5 n6 n7 x1 x2 tbl chi2stat pval

% piefig=figure; piefig.Color=[1 1 1];
% ax=subplot(1,2,1);
% p1=pie([numCorC(1),numCorC(2),numCorC(3)]);
% title(ax,'Control');ax.FontSize=20;
% hold on
% ax=subplot(1,2,2);
% p2=pie([numCorP(1),numCorP(2),numCorP(3)]);
% title(ax,'PAE');ax.FontSize=20;
% print(piefig,'-bestfit', '-dpdf', '-r600',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/Learning&Memory/',filesep,'PIEcorrs_fig.pdf'])


% RATE OVERLAP
% AllStats=ScatterBox([ratediff1C;ratediff2C],[ratediff1P;ratediff2P],{'Sacc','PAE'},{'Rate Overlap'},2)
% AllStats=CDFplots([ratediff1C;ratediff2C],[ratediff1P;ratediff2P],{'Sacc','PAE'},{'Rate Overlap'},2)
% print(figure(1),'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/Gordon Conference/measure_figs',filesep,'Rate Overlap_fig.pdf'])

% print(figure(1),'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub',filesep,'rateoverlap_fig.pdf'])
end
%%  THETA MODULATION
[thetaindexoutC,peakoutC,coroutC]=collectautocorrs([ids1C(~ismember(nids1C,nids2C),:);ids2C],data);
[thetaindexoutP,peakoutP,coroutP]=collectautocorrs([ids1P(~ismember(nids1P,nids2P),:);ids2P],data);

AllStats=CDFplots(thetaindexoutC,thetaindexoutP,{'Sacc','PAE'},{'Theta Index'},2)

[s,I]=sort(thetaindexoutC);
sortautocorC=coroutC(I,:);

[s,I]=sort(thetaindexoutP);
sortautocorP=coroutP(I,:);

% tempID=[ids1P(~ismember(nids1P,nids2P),:);ids2P];
% tempID(I,:);
% 
% 
% [thetaindex,peak,cor,lag] = thetamodulation(data.RH11.S20160520111231.Spikes{20, 1}/1000000 );

%% PICK GOOD CELLS FOR DISPLAY
pickcells=1;
id=ids1C;
if pickcells==1
    for i=1:size(id,1)
        % EXTRACT FRAME AND SPIKE DATA TO PLOT TRACK
        frames=data.(id{i,1}).(id{i,2}).frames;
        % get spikes
        cellnum=regexp(id(i,3),'\d*','Match');
        spkts=data.(id{i,1}).(id{i,2}).Spikes{str2double(cellnum{1})};
        % get events
        events=data.(id{i,1}).(id{i,2}).events;
        
        % find track length
        date=char(id{i,2});
        date=[date(2:5),'-',date(6:7),'-',date(8:9),'_',date(10:11),'-',date(12:13),'-',date(14:15)];
        track_length=TrackLength(['/Volumes/Ryan_4TB/Place_Cell_Data/RawPAE_PlaceCell',filesep,id{i,1},filesep,date]);
        
        [right,left,stability,time2stableC,firststablelap,meanstability,lap_stability]=rebuildFramematrix(frames,spkts,events,track_length);
        clear date events spkts cellnum frames
        
        fig=figure;
        fig.Color=[1 1 1];
        fig.OuterPosition=[564 6 877 799];
        
        subplot(2,2,1)
        plot(right.dataspks(:,2),right.dataspks(:,1),'.k');hold on
        scatter(right.dataspks(right.dataspks(:,6)==1,2),right.dataspks(right.dataspks(:,6)==1,1),'filled','r')
        ylabel('Laps')
        xlabel('Distance on track (cm)')
        set(gca,'FontWeight','bold','FontSize',18,'LineWidth',3,'XTick',...
            linspace(min(right.dataspks(:,2)),max(right.dataspks(:,2)),3),...
            'XTickLabel',[0 60 120],'box','off','YTickLabel',[],'TickLength',[0 0])
        xlim([min(right.dataspks(:,2)),max(right.dataspks(:,2))])
        ylim([min(right.dataspks(:,1)),max(right.dataspks(:,1))])
        
        
        subplot(2,2,2)
        plot(left.dataspks(:,2),left.dataspks(:,1),'.k');hold on
        scatter(left.dataspks(left.dataspks(:,6)==1,2),left.dataspks(left.dataspks(:,6)==1,1),'filled','r')
        ylabel('Laps')
        xlabel('Distance on track (cm)')
        set(gca,'FontWeight','bold','FontSize',18,'LineWidth',3,'XTick',...
            linspace(min(left.dataspks(:,2)),max(left.dataspks(:,2)),3),...
            'XTickLabel',[0 60 120],'box','off','YTickLabel',[],'TickLength',[0 0])
        xlim([min(left.dataspks(:,2)),max(left.dataspks(:,2))])
        ylim([min(left.dataspks(:,1)),max(left.dataspks(:,1))])
        
        subplot(2,2,3)
        h=rainbowplot(right.SmoothRateMap,0);
        set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3,...
            'XTick',[1 length(right.SmoothRateMap)/2 length(right.SmoothRateMap)],...
            'XTickLabel',[0 60 120],'TickLength',[0 0])
        ylabel('Firing Rate (hz)')
        xlabel('Distance on track (cm)')
        ylim([0 max([right.SmoothRateMap,left.SmoothRateMap])])
        hold off
        
        subplot(2,2,4)
        h=rainbowplot(left.SmoothRateMap,0);
        set(gca,'box','off','FontWeight','bold','FontSize',18,'LineWidth',3,...
            'XTick',[1 length(left.SmoothRateMap)/2 length(left.SmoothRateMap)],...
            'XTickLabel',[0 60 120],'TickLength',[0 0])
        ylabel('Firing Rate (hz)')
        xlabel('Distance on track (cm)')
        ylim([0 max([right.SmoothRateMap,left.SmoothRateMap])])
        hold off
        %                 print(fig,'-dpng', '-r600',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/Learning&Memory/ControlExamples',filesep,id{i,:},'.png'])
        
        close all
    end
end


%% SAMPLING BIAS QUESTION
samplingbias=0;
if samplingbias==1
    for ii=1:10000
        disp([num2str(ii),' of 10000'])
        [y,idx]= datasample(paemap,232);
        [s,I]=sort(idx);
        tempmat=y(I,:);
        TL=tempmat(1:116,1:60);
        TR=tempmat(1:116,61:end);
        BL=tempmat(117:end,1:60);
        BR=tempmat(117:end,61:end);
        tempcorrmat=[[corr(TL,TL),corr(TL,TR)];[corr(BL,BR),corr(BR,BR)]];
        means=[];
        sems=[];
        for i=0:59
            tempvec=[diag(tempcorrmat(61:end,1:60),i);diag(tempcorrmat(61:end,1:60),-i);diag(tempcorrmat(1:60,61:end),i);diag(tempcorrmat(1:60,61:end),-i)];
            means(i+1)=mean(tempvec);
            sems(i+1)=std(tempvec)/sqrt(length(tempvec));
        end
        storemean(ii,:)=means;
        storsem(ii,:)=sems;
    end
    mean(storemean,1);
end
%%

% plot right and left in one matrix
% control
popvecfig=figure;
popvecfig.Color=[1 1 1];

controlmap=[mapright_1C,mapleft_1C;mapright_2C,mapleft_2C];
subplot(1,2,1)
p1=imagesc(controlmap);colormap jet
set(gca,'Visible','off','box','off')

paemap=[mapright_1P,mapleft_1P;mapright_2P,mapleft_2P];
subplot(1,2,2)
p2=imagesc(paemap);colormap jet
set(gca,'Visible','off','box','off')
% print(popvecfig,'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub',filesep,'popvecfig.pdf'])


popveccor_fig=figure;
popveccor_fig.Color=[1 1 1];
controlauto=[[corr(mapright_1C,mapright_1C);corr(mapright_1C,mapleft_1C)],[corr(mapright_2C,mapleft_2C);corr(mapleft_2C,mapleft_2C)]];
paeauto=[[corr(mapright_1P,mapright_1P);corr(mapright_1P,mapleft_1P)],[corr(mapright_2P,mapleft_2P);corr(mapleft_2P,mapleft_2P)]];
subplot(1,2,1)
imagesc(controlauto);colormap jet;set(gca,'Visible','off','box','off');axis image
subplot(1,2,2)
imagesc(paeauto);colormap jet;set(gca,'Visible','off','box','off');axis image
% print(popveccor_fig,'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub',filesep,'popveccor_fig.pdf'])



%% decor distance diag method
decordistance_fig=figure;
decordistance_fig.Color=[1 1 1];
means=[];
sems=[];
for i=0:59
    tempvec=[diag(controlauto(61:end,1:60),i);diag(controlauto(61:end,1:60),-i);diag(controlauto(1:60,61:end),i);diag(controlauto(1:60,61:end),-i)];
    means(i+1)=mean(tempvec);
    sems(i+1)=std(tempvec)/sqrt(length(tempvec));
end
decorc=means;
varargout=shadedErrorBar(linspace(1,60,length(means)),means,sems,'-k',1);hold on

for i=0:59
    tempvec=[diag(paeauto(61:end,1:60),i);diag(paeauto(61:end,1:60),-i);diag(paeauto(1:60,61:end),i);diag(paeauto(1:60,61:end),-i)];
    means(i+1)=mean(tempvec);
    sems(i+1)=std(tempvec)/sqrt(length(tempvec));
end
decorp=means;
varargout=shadedErrorBar(linspace(1,60,length(means)),means,sems,'-r',1);
set(gca,'box','off','FontWeight','bold','FontSize',25,'LineWidth',3)
xlabel(['Distance (cm)'])
ylabel(['Correlation'])

AllStats=ScatterBox(decorc',decorp',{'Sacc','PAE'},{'Decor Distance'},2)

% print(decordistance_fig,'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub',filesep,'decordistance_fig.pdf'])



% decor distance of entire corr vector
%     fulldecordistance_fig=figure;
%     fulldecordistance_fig.Color=[1 1 1];
%     for i=1:120
%         tempvec=[diag(controlauto,i);diag(controlauto,-i)];
%         means(i)=mean(tempvec);
%         sems(i)=std(tempvec)/sqrt(length(tempvec));
%     end
%     varargout=shadedErrorBar(linspace(1,60,length(means)),means,sems,'-k',1);hold on
%
%     for i=1:120
%         tempvec=[diag(paeauto,i);diag(paeauto,-i)];
%         means(i)=mean(tempvec);
%         sems(i)=std(tempvec)/sqrt(length(tempvec));
%     end
%     varargout=shadedErrorBar(linspace(1,60,length(means)),means,sems,'-r',1);
%     set(gca,'box','off','FontWeight','bold','FontSize',25,'LineWidth',3)
%     xlabel(['Distance (cm)'])
%     ylabel(['Correlation'])
% print(fulldecordistance_fig,'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub',filesep,'diagcorr_fig.pdf'])


%% DOWN DIAG CORR VALUES
% DIAG 0
centerdiag_fig=figure;
centerdiag_fig.Color=[1 1 1];
controldiag=[diag(controlauto(61:end,1:60),0),diag(controlauto(1:60,61:end),0)];
varargout=shadedErrorBar(linspace(1,60,length(controldiag)),mean(controldiag,2),std(controldiag,[],2)/sqrt(2),'-k',1);hold on

paediag=[diag(paeauto(61:end,1:60),0),diag(paeauto(1:60,61:end),0)];
varargout=shadedErrorBar(linspace(1,60,length(paediag)),mean(paediag,2),std(paediag,[],2)/sqrt(2),'-r',1);hold on
set(gca,'box','off','FontWeight','bold','FontSize',25,'LineWidth',3)
xlabel(['Distance (cm)'])
ylabel(['Correlation'])

AllStats=ScatterBox(controldiag,paediag,{'Sacc','PAE'},{'DiagCorr'},2)
% print(centerdiag_fig,'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub',filesep,'centerdiag_fig.pdf'])


%% CUMULATIVE FREQ OF ROW BY ROW CORRELATIONS
figure;
% subplot(1,2,1)
[f,x] = ecdf([RowCorr_1C;RowCorr_2C]);plot(x,f,'DisplayName','Control','LineWidth',4,'Color','k');hold on
[f,x] = ecdf([RowCorr_1P;RowCorr_2P]);plot(x,f,'DisplayName','PAE','LineWidth',4,'Color','r')
box off;legend('show','Location','best')
xlabel('Directional Correlation');ylabel('Cumulative Frequency')

AllStats=CDFplots([RowCorr_1C;RowCorr_2C],[RowCorr_1P;RowCorr_2P],{'Sacc','PAE'},{'Correlation'},3)
% print(figure(1),'-bestfit', '-dpdf', '-r600',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub',filesep,'SpatialCorrelations_fig.pdf'])


% \/ LOCAL FUNCTIONS \/
% /////////////////////////////////////////////////////////////////////////////////////////////////////////
%% CHECK AND STORE NORMALIZED PLACE CELLS
function [mapright,mapleft,RowCorr,ratediff,unnormright,unnormleft,ids]=extract(rats,data,direc)
mapright=[];
mapleft=[];
RowCorr=[];
ratediff=[];
unnormright=[];
unnormleft=[];
ids=[];
for r=1:length(rats)
    sessions=fieldnames(data.(rats{r}));
    for s=1:length(sessions)
        cells=fieldnames(data.(rats{r}).(sessions{s}).ratemap);
        %
        %         cells=cells(data.(rats{r}).(sessions{s}).measures(:,1,direc)>0.25 &... % info content >.25
        %             data.(rats{r}).(sessions{s}).measures(:,5,direc)>0.3 &... % >.0 hz average rate
        %             data.(rats{r}).(sessions{s}).measures(:,8,direc)>50 &... % >50 spikes
        %             data.(rats{r}).(sessions{s}).measures(:,38,direc)<.05 &... % <0.5% spikes <2ms in autocorr
        %             data.(rats{r}).(sessions{s}).measures(:,39,direc)<0.5 &... % <50% cut off by threshold
        %             data.(rats{r}).(sessions{s}).measures(:,5,direc)<10 &... % <10hz avg firing
        %             data.(rats{r}).(sessions{s}).measures(:,48,direc)>0.2 &...% >.2ms peak to valley duration
        %             data.(rats{r}).(sessions{s}).measures(:,54,direc)>5 &...% above 5 laps
        %             data.(rats{r}).(sessions{s}).measures(:,53,direc)>.02); % temporal stability >.02
        
        cells=cells(data.(rats{r}).(sessions{s}).measures(:,1,direc)>0.25 &... % info content >.25
            data.(rats{r}).(sessions{s}).measures(:,5,direc)>0.3 &... % >.3 hz average rate
            data.(rats{r}).(sessions{s}).measures(:,4,direc)>2 &...
            data.(rats{r}).(sessions{s}).measures(:,8,direc)>50 &... % >50 spikes
            data.(rats{r}).(sessions{s}).measures(:,38,direc)<.05 &... % <0.5% spikes <2ms in autocorr
            data.(rats{r}).(sessions{s}).measures(:,39,direc)<0.5 &... % <50% cut off by threshold
            data.(rats{r}).(sessions{s}).measures(:,54,direc)>5 &...% above 5 laps
            data.(rats{r}).(sessions{s}).measures(:,53,direc)>.02); % temporal stability >.02
        
        %                 cells=cells(data.(rats{r}).(sessions{s}).measures(:,1,direc)>1 &...
        %                     data.(rats{r}).(sessions{s}).measures(:,8,direc)>200 &...
        %                     data.(rats{r}).(sessions{s}).measures(:,31,direc)>data.(rats{r}).(sessions{s}).measures(:,8,direc)*.2 &... % >20% of spikes in firing field
        %                     data.(rats{r}).(sessions{s}).measures(:,24,1)>0,:,:);
        
        %         cells=cells(data.(rats{r}).(sessions{s}).measures(:,1,direc)>.6 &...
        %             data.(rats{r}).(sessions{s}).measures(:,8,direc)>100 &...
        %             data.(rats{r}).(sessions{s}).measures(:,21,direc)>10,:,:);
        
        for c=1:length(cells)
            if ~isfield(data.(rats{r}).(sessions{s}).ratemap.(cells{c}),'session1')
                continue
            end
            ids=[ids;{(rats{r}),(sessions{s}),(cells{c})}];
            
            % resize to 60 bins
            map1=imresize(data.(rats{r}).(sessions{s}).ratemap.(cells{c}).session1,[1,60]);
            map2=imresize(data.(rats{r}).(sessions{s}).ratemap.(cells{c}).session2,[1,60]);
            
            % norm rate maps together
            normed=rescale([map1,map2],0,1);
            map1=normed(1:60);
            map2=normed(61:end);
            
            % collect ratemaps
            mapright=[mapright;map1];
            
            mapleft=[mapleft;map2];
            
            % collect correlations
            RowCorr=[RowCorr;corr2(map1,map2)];
            
            % collect rate overlap
            if direc==1
                ratediff=[ratediff;(max(map1)-max(map2))];
            elseif direc==2
                ratediff=[ratediff;(max(map2)-max(map1))];
            end
            
            % collect un-normalized ratemaps
            unnormright=[unnormright;imresize(data.(rats{r}).(sessions{s}).ratemap.(cells{c}).session1,[1,60])];
            unnormleft=[unnormleft;imresize(data.(rats{r}).(sessions{s}).ratemap.(cells{c}).session2,[1,60])];
        end
    end
end
end

% ARRANGE BINS RIGHT
function [mat1,mat2]=arrange(mat1,mat2)
[~,I]=max(mat1,[],2);
[~,I2]=sort(I);
mat1=mat1(I2,:);
mat2=mat2(I2,:);
end

function cor=flipcorr(mat1,mat2)
cor=[];
for i=1:length(mat1)
    cor=[cor;corr2(mat1(i,:),flip(mat2(i,:),2))];
end
end

function matfordis=autocrosscorr(matr,matl)
autoR=corr(matr,matr);
R=corr(matr,matl);
matfordis=[autoR,R;R,autoR];
end



function [right,left,stability,time2stable,firststablelap,meanstability,lap_stability]=rebuildFramematrix(frames,SpikeFile,events,track_length)
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


% % VELO FILTER BASED ON INDEX CREATED ABOVE
% data_video_nospk(logical(in),:)=[];
% data_video_nospk(:,6)=[];
% data_video_spk(data_video_spk(:,6)==1,:)=[];
% data_video_spk(:,6)=[];

% SPLIT RIGHT AND LEFT RUNS
[right,left,DirectionalityIndex,Displacement,nlaps]=RightVsLeft(data_video_nospk,data_video_spk,track_length,30);


[lapsR,mapsR,corsR]=splitlapagain(right.dataspks,find(diff(right.dataspks(:,1))>1e7),track_length,right.SmoothRateMap);
[lapsL,mapsL,corsL]=splitlapagain(left.dataspks,find(diff(left.dataspks(:,1))>1e7),track_length,left.SmoothRateMap);

corrz=[];
for m=1:size(mapsR,1)
    for mm=1:size(mapsR,1)
        corrz=[corrz;corr2(mapsR(m,:),mapsR(mm,:))];
    end
end
corrz(1:size(mapsR,1)+1:length(corrz))=[];
lap_perm_stabilityR=nanmean(corrz);


corrz=[];
for m=1:size(mapsL,1)
    for mm=1:size(mapsL,1)
        corrz=[corrz;corr2(mapsL(m,:),mapsL(mm,:))];
    end
end
corrz(1:size(mapsL,1)+1:length(corrz))=[];
lap_perm_stabilityL=nanmean(corrz);

lap_stability=nanmean([lap_perm_stabilityR,lap_perm_stabilityL]);

[ICright]=infocontent(right.SmoothRateMap,right.occ,right.nBinsx,right.nBinsy);
[ICleft]=infocontent(left.SmoothRateMap,left.occ,left.nBinsx,left.nBinsy);

if ICright>.25 && max(right.SmoothRateMap)>2 
    rightplacecell=1;
else
    rightplacecell=0;
end
if ICleft>.25 && max(left.SmoothRateMap)>2
    leftplacecell=1;
else
    leftplacecell=0;
end


if rightplacecell==1 && leftplacecell==1
    lap_stability=nanmean([lap_perm_stabilityR,lap_perm_stabilityL]);
    meanstability=nanmean([corsR;corsL]);
elseif rightplacecell==1 && leftplacecell==0
    lap_stability=lap_perm_stabilityR;
    meanstability=nanmean(corsR);
elseif rightplacecell==0 && leftplacecell==1
    lap_stability=lap_perm_stabilityL;
    meanstability=nanmean(corsL);
elseif rightplacecell==0 && leftplacecell==0
    lap_stability=nanmean([lap_perm_stabilityR,lap_perm_stabilityL]);
    meanstability=nanmean([corsR;corsL]);
end

if ~exist('lap_stability','var') || ~exist('meanstability','var') 
   test=1 
end

% meanstability=nanmean([corsR;corsL]);


stability=[corr2([1:length(corsL)]',corsL),corr2([1:length(corsR)]',corsR)];

% find time to stability
stablelap=find((corsR>.2));
if isempty(stablelap)
    time2stable(1,2)=NaN;
    firststablelap(1,2)=NaN;
else
    neartime2stable=find(data_video_nospk(:,1)>=lapsR{stablelap(1)}(1,1));
    time2stable(1,2)=neartime2stable(1)/30;
    firststablelap(1,2)=stablelap(1);
end

stablelap=find((corsL>.2));
if isempty(stablelap)
    time2stable(1,1)=NaN;
    firststablelap(1,1)=NaN;
else
    neartime2stable=find(data_video_nospk(:,1)>=lapsL{stablelap(1)}(1,1));
    time2stable(1,1)=neartime2stable(1)/30;
    firststablelap(1,1)=stablelap(1);
end

% plot option
plotz=0;
if plotz==1
    stabFig=figure;stabFig.Color=[1 1 1];
    subplot(1,3,1)
    plot(right.dataspks(:,2),right.dataspks(:,1),'.k');hold on
    scatter(right.dataspks(right.dataspks(:,6)==1,2),right.dataspks(right.dataspks(:,6)==1,1),'filled','r')
    ylabel('Laps')
    xlabel('Distance on track (cm)')
    set(gca,'FontWeight','bold','FontSize',18,'LineWidth',3,'XTick',...
        linspace(min(right.dataspks(:,2)),max(right.dataspks(:,2)),3),...
        'XTickLabel',[0 60 120],'box','off','YTickLabel',[],'TickLength',[0 0])
    xlim([min(right.dataspks(:,2)),max(right.dataspks(:,2))])
    ylim([min(right.dataspks(:,1)),max(right.dataspks(:,1))])
    
    subplot(1,3,2)
    tempmaps=mapsR;
    for i=1:size(tempmaps,1)
        rmp=plot(rescale(tempmaps(i,:),0,10000000)+repmat(mean(lapsR{i}(:,1)),1,size(tempmaps,2)),'k');hold on
        rmp.LineWidth=2;
    end
    xlabel('Distance on track (cm)')
    set(gca,'FontWeight','bold','FontSize',18,'LineWidth',3,'XTick',...
        [1 length(right.SmoothRateMap)/2 length(right.SmoothRateMap)],...
        'XTickLabel',[0 60 120],'box','off','YTickLabel',[],'TickLength',[0 0])
    xlim([1,size(tempmaps,2)])
    ylim([lapsR{1}(1,1),max(rescale(tempmaps(i,:),0,10000000)+repmat(mean(lapsR{i}(:,1)),1,size(tempmaps,2)))])
    ylabel('Laps')
    
    subplot(1,3,3);
    for i=1:size(tempmaps,1)
        meantime(i,1)=mean(lapsR{i}(:,1));
    end
    p = polyfit(corsR,meantime,1);
    f = polyval(p,corsR);
    scatter(corsR,meantime,'Filled','r');hold on
    lsp=plot(corsR,f,'-');
    lsp.Color='k';
    lsp.LineWidth=2;
    set(gca,'FontWeight','bold','FontSize',18,'LineWidth',3,'box','off','YTickLabel',[],'TickLength',[0 0])
    ylim([lapsR{1}(1,1),meantime(end)])
    xlabel('Correlation (r)')
    ylabel('Laps')
    %     print(stabFig,'-dpng', '-r600',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/Learning&Memory/',filesep,'stabFigExample.png'])
end
end

function [lap,SmoothRateMap,cors]=splitlapagain(tempmat,idx,track_length,map)
i=1;
while true
    if i>length(idx)
        lap{i,1}=tempmat(idx(i-1)+1:end,:);
        [SmoothRateMap(i,:)]=binlaps(vertcat(lap{:,:}),lap{i,1}((lap{i,1}(:,6)==0),:),...
            30,lap{i,1}((lap{i,1}(:,6)==1),:),track_length);
        for m=1:size(SmoothRateMap,1)
            cors(m,1)=corr2(map,(SmoothRateMap(m,:)));
        end
        cors(isnan(cors))=0;
        break
    end
    if i==1
        lap{i,1}=tempmat(1:idx(i),:);
        [SmoothRateMap(i,:)]=binlaps(vertcat(lap{:,:}),lap{i,1}((lap{i,1}(:,6)==0),:),...
            30,lap{i,1}((lap{i,1}(:,6)==1),:),track_length);
    else
        lap{i,1}=tempmat(idx(i-1)+1:idx(i),:);
        [SmoothRateMap(i,:)]=binlaps(vertcat(lap{:,:}),lap{i,1}((lap{i,1}(:,6)==0),:),...
            30,lap{i,1}((lap{i,1}(:,6)==1),:),track_length);
    end
    i=i+1;
end
end

function [stability,time2stable,firststablelap,meanstability,lap_stability]=betweenlapstability(id,data)
for i=1:size(id,1)
    % EXTRACT FRAME AND SPIKE DATA TO PLOT TRACK
    frames=data.(id{i,1}).(id{i,2}).frames;
    % get spikes
    cellnum=regexp(id(i,3),'\d*','Match');
    spkts=data.(id{i,1}).(id{i,2}).Spikes{str2double(cellnum{1})};
    % get events
    events=data.(id{i,1}).(id{i,2}).events;
    
    % find track length
    date=char(id{i,2});
    date=[date(2:5),'-',date(6:7),'-',date(8:9),'_',date(10:11),'-',date(12:13),'-',date(14:15)];
    track_length=TrackLength(['/Volumes/Ryan_4TB/Place_Cell_Data/RawPAE_PlaceCell',filesep,id{i,1},filesep,date]);
    
    [~,~,stability(i,:),time2stable(i,:),firststablelap(i,:),meanstability(i,:),lap_stability(i,:)]=rebuildFramematrix(frames,spkts,events,track_length);
    clear date events spkts cellnum frames
end
end

function [SmoothRateMap]=binlaps(occ_overall,occMatrix,sampleRate,spks_VEL,track_length)

nBinsx = round(track_length/2); nBinsy = 1;

if isempty(occMatrix)
    SmoothRateMap=zeros(nBinsy,nBinsx);
    return
end

MinY = min(occ_overall(:,3));
MaxY = max(occ_overall(:,3));
MinX = min(occ_overall(:,2));
MaxX = max(occ_overall(:,2));
edges{1} = linspace(MinY, MaxY, nBinsy+1);
edges{2} = linspace(MinX, MaxX, nBinsx+1);

occMatrix = [occMatrix(:,3),occMatrix(:,2)];
Omatrix = hist3(occMatrix,'Edges',edges);
Omatrix(2,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/sampleRate;

% bin spike data
Smatrix = hist3([spks_VEL(:,3), spks_VEL(:,2)],'Edges',edges);
Smatrix(2,:) = [];
Smatrix(:,end) = [];

% divide binned spikes by occupancy to get rate maps and store data
% (can use occ+eps instead of removing 0's)
FilledRateMatrix = Smatrix./occ;
FilledRateMatrix(isnan(FilledRateMatrix)) = 0;
FilledRateMatrix(isinf(FilledRateMatrix)) = 0;
% FilledRateMatrix(occ<0.150)=0;

% SMOOTH
filtWidth = [1,5]; filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
end

function [fieldoverlap,Displacement,rateoverlap,distalcode,positioncode]=OVERLAP(id,data)
track_length=120;
for i=1:size(id,1)
    % EXTRACT RATEMAP
    cellnum=regexp(id(i,3),'\d*','Match');
    ratemap1=data.(id{i,1}).(id{i,2}).ratemap.(char(strcat('Cell',cellnum{1}))).session1;
    ratemap2=data.(id{i,1}).(id{i,2}).ratemap.(char(strcat('Cell',cellnum{1}))).session2;
    
    if length(ratemap1)<60
        track_length=90;
    end
    
    % FIELD OVERLAP
    [~,peak1]=max(ratemap1);
    [~,peak2]=max(ratemap2);
    if peak1==peak2
        fieldoverlap(i,1)=1;
    else
        [field1]=findfield(ratemap1);
        [field2]=findfield(ratemap2);
        fieldoverlap(i,1)=length(intersect(find(field1),find(field2)))/sum([field1+field2]>0);
    end
    
    % DISPLACEMENT
    [~,i1]=max(ratemap1);[~,i2]=max(ratemap2);
    Displacement(i,1)=abs(i1-i2)*(track_length/length(ratemap1));
    
    % RATE OVERLAP
    normratemaps=rescale([ratemap1,ratemap2],0,1);
    norm1=normratemaps(1:length(ratemap1));
    norm2=normratemaps(length(ratemap1)+1:end);
    rateoverlap(i,1)=abs(max(norm1)-max(norm2));
    
    %     DirectionalityIndex(i,1)=abs(sum(ratemap1-ratemap2)/sum(ratemap1+ratemap2));
    
    % DISTAL CODING
    tempmap=flip(ratemap1,2);
    [m,I]=max(tempmap);
    if I==i2 || I+1==i2 || I-1==i2
        distalcode(i,1)=1;
    else
        distalcode(i,1)=0;
    end
    
    % POSITION
    if i1==i2 || i1==i2+1 || i1==i2-1
        positioncode(i,1)=1;
    else
        positioncode(i,1)=0;
    end
    %     figure;plot(ratemap1);hold on;plot(ratemap2)
    %     close all
end
end

function [field]=findfield(Ratemap)
[M,I]=max(Ratemap); Thres=M*0.20;
field=zeros(1,length(Ratemap)); field(I)=1;

% FORWARD
for i=1:length(Ratemap)-I
    if Ratemap(I)~=length(Ratemap)
        if Ratemap(I+i)>Thres
            field(I+i)=1;
        elseif Ratemap(I+i)<Thres
            break
        end
    end
end

% BACKWARD
for i=1:I-1
    if Ratemap(I)~=1
        if Ratemap(I-i)>Thres
            field(I-i)=1;
        elseif Ratemap(I-i)<Thres
            break
        end
    end
end
end

function [idxmat]=labeldirectionalcodes(positioncodeC,distalcodeC,rateoverlapC)
% justfieldC,justrateC,justdistalC,justpositionC,ratebyfieldC,ratebydistalC,unidentified
% bidirectional:
%   1: change positions
%   2: can't code for distance
%   3: can't change in rate over 20%
justfieldC=sum(positioncodeC==0 & distalcodeC==0 & rateoverlapC<=.20);
idxmat(:,1)=[positioncodeC==0 & distalcodeC==0 & rateoverlapC<=.20]';
% Rate code
%   1: Must have >20% rate change
%   2: must be in same position
justrateC=sum(rateoverlapC>.20 & positioncodeC==1 & distalcodeC==0);
idxmat(:,2)=[rateoverlapC>.20 & positioncodeC==1 & distalcodeC==0]';
% distal code
%   1: must mirror itself
%   2: must not rate change
justdistalC=sum(distalcodeC==1 & rateoverlapC<=.20 & positioncodeC==0);
idxmat(:,3)=[distalcodeC==1 & rateoverlapC<=.20 & positioncodeC==0]';
% position code
%   1: must be in same position
%   2: no rate change
justpositionC=sum(positioncodeC==1 & rateoverlapC<=.20 & distalcodeC==0);
idxmat(:,4)=[positioncodeC==1 & rateoverlapC<=.20 & distalcodeC==0]';
% bidirectional by Rate
%   1: change positions
%   2: can't code for distance
%   3: Change in rate over 20%
ratebyfieldC=sum(positioncodeC==0 & distalcodeC==0 & rateoverlapC>.20);
idxmat(:,5)=[positioncodeC==0 & distalcodeC==0 & rateoverlapC>.20]';
% distal by rate
%   1: must mirror itself
%   2: rate change
ratebydistalC=sum(distalcodeC==1 & rateoverlapC>.20 & positioncodeC==0);
idxmat(:,6)=[distalcodeC==1 & rateoverlapC>.20 & positioncodeC==0]';

idxmat((sum(idxmat,2)==0),7)=1;

% sum([justfieldC,justrateC,justdistalC,justpositionC,ratebyfieldC,ratebydistalC]);
end

function [thetaindexout,peakout,corout]=collectautocorrs(id,data)
for i=1:length(id)
%     cellnum=regexp(id(i,3),'\d*','Match');
%     autocor(i,:)=data.(id{i,1}).(id{i,2}).thetaautocorr.(char(strcat('Cell',cellnum{1}))).session1;
    
    % EXTRACT FRAME AND SPIKE DATA TO PLOT TRACK
    frames=data.(id{i,1}).(id{i,2}).frames;
    % get spikes
    cellnum=regexp(id(i,3),'\d*','Match');
    spkts=data.(id{i,1}).(id{i,2}).Spikes{str2double(cellnum{1})};
    % get events
    events=data.(id{i,1}).(id{i,2}).events;
    
    % find track length
    date=char(id{i,2});
    date=[date(2:5),'-',date(6:7),'-',date(8:9),'_',date(10:11),'-',date(12:13),'-',date(14:15)];
    track_length=TrackLength(['/Volumes/Ryan_4TB/Place_Cell_Data/RawPAE_PlaceCell',filesep,id{i,1},filesep,date]);
    
    [right,left,stability,time2stable,firststablelap,meanstability]=rebuildFramematrix(frames,spkts,events,track_length);
    frames=[right.datavid;left.datavid];
    ts=sort(frames(:,1));
    ts=unique(ts);
    
    spk=[right.dataspks;left.dataspks];
    spk=spk(spk(:,6)==1,1);
    spk=sort(spk(:,1));
    
    sec=linspace(0,(length(ts)/30),length(ts));
 
    secspks=interp1(ts',sec',spk);
    
    [thetaindex,peak,cor,lag] = thetamodulation(secspks);
%     figure;plot(lag,cor)
%     
%         grade=data.(id{i,1}).(id{i,2}).measures(str2double(cellnum{1}),24,1)
%         isi=data.(id{i,1}).(id{i,2}).measures(str2double(cellnum{1}),38,1)

thetaindexout(i,1)=thetaindex;
peakout(i,1)=peak;
corout(i,:)=cor;
end
end

function [InformationContent]=infocontent(SmoothRateMap,occ,nBinsx,nBinsy)
% calculate information content from nsma
rY=reshape(SmoothRateMap,nBinsx*nBinsy,1);rY(isnan(rY))=0;rY(isinf(rY))=0;
occRSHP=reshape(occ,nBinsx*nBinsy,1);occSUM=sum(occRSHP);pX=occRSHP./occSUM;
[nBins,nCells]=size(rY);relR=rY./kron(ones(nBins,1),pX'*rY);
log_relR=log2(relR);log_relR(isinf(log_relR))=0;
InformationContent=sum(kron(pX,ones(1,nCells)).*relR.*log_relR);
end

%                                                      ___
%                                                   ,o88888
%                                                ,o8888888'
%                          ,:o:o:oooo.        ,8O88Pd8888"
%                      ,.::.::o:ooooOoOoO. ,oO8O8Pd888'"
%                    ,.:.::o:ooOoOoOO8O8OOo.8OOPd8O8O"
%                   , ..:.::o:ooOoOOOO8OOOOo.FdO8O8"
%                  , ..:.::o:ooOoOO8O888O8O,COCOO"
%                 , . ..:.::o:ooOoOOOO8OOOOCOCO"
%                  . ..:.::o:ooOoOoOO8O8OCCCC"o
%                     . ..:.::o:ooooOoCoCCC"o:o
%                     . ..:.::o:o:,cooooCo"oo:o:
%                  `   . . ..:.:cocoooo"'o:o:::'
%                  .`   . ..::ccccoc"'o:o:o:::'
%                 :.:.    ,c:cccc"':.:.:.:.:.'
%               ..:.:"'`::::c:"'..:.:.:.:.:.'
%             ...:.'.:.::::"'    . . . . .'
%            .. . ....:."' `   .  . . ''
%          . . . ...."'
%          .. . ."'
%         .


% %% CREATE GIF
% gifmaker=0;
% if gifmaker==1
%     % CREATE SESSION GIF
%     ids1C(56,:)
%     % get video
%     frames=data.(ids1C{56,1}).(ids1C{56,2}).frames;
%     % get spikes
%     cellnum=regexp(ids1C(56,3),'\d*','Match');
%     spkts=data.(ids1C{56,1}).(ids1C{56,2}).Spikes{str2double(cellnum{1})};
%     X=interp1(frames(:,1),frames(:,2),spkts);
%     Y=interp1(frames(:,1),frames(:,3),spkts);
%
%     data_video_spk=sortrows([[spkts,X,Y ones(size(X,1),1)];[frames,zeros(length(frames),1)]],1);
%
%     gif_fig=figure;
%     gif_fig.Color=[1 1 1];
%     gif_fig.OuterPosition=[80 463 1268 187];
%     xlim([min(data_video_spk(:,2)) max(data_video_spk(:,2))])
%     ylim([min(data_video_spk(:,3)) max(data_video_spk(:,3))])
%     linearbox=[[min(data_video_spk(:,2)),min(data_video_spk(:,3))];...
%         [min(data_video_spk(:,2)),max(data_video_spk(:,3))];...
%         [max(data_video_spk(:,2)),max(data_video_spk(:,3))];...
%         max(data_video_spk(:,2)),min(data_video_spk(:,3));...
%         [min(data_video_spk(:,2)),min(data_video_spk(:,3))]];
%     plot(linearbox(:,1),linearbox(:,2),'LineWidth',3,'color','k')
%     box off
%     axis off
%     axis tight manual
%     hold on
%
%     % downsampled method
%     plot(data_video_spk(:,2),data_video_spk(:,3),'w');
%     hold on
%     p = plot(data_video_spk(1,2),data_video_spk(1,3),'square','MarkerFaceColor','k','MarkerEdgeColor','k','Linewidth',40);
%     % hold off
%     axis manual
%     SXY=data_video_spk;
%     SXY(SXY(:,4)~=1,[2:3])=NaN;
%     nframe=1;
%     for k = 2:60:length(data_video_spk(:,2))
%         p.XData = data_video_spk(k,2);
%         p.YData = data_video_spk(k,3);
%         drawnow
%         if k-60<0;fac=2;else;fac=k;fac=fac-60;end
%         scatter(SXY([fac:k],2),SXY([fac:k],3),50,'filled','r')
%         uistack(p,'top');
%         print(gif_fig,'-djpeg',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub/gif2',filesep,['gif_fig',num2str(nframe),'.jpeg']])
%         nframe=nframe+1;
%     end
% end

