% intra_trial_stability
% REVIEWER 2
% Comment 2: If the tilted mice have path integration problems, this would
% be expected to show up as a reduction in intra-trial stability:
% the authors should report a temporal stability measure, for example
% stability across time-slices of a trial, to look for this.

% REVIEWER 3
% Comment 5: A measure of stability should be added in all sessions -
% check the correlation between the 1st and 2nd half of each session and
% see whether there is a difference between controls and titled mice.
% This should help the reader comprehend the issue of stability better in this paper.

clear;clc;close all
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewOCC_Map_workspace2.mat')
FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';
 
resultC=analyze(dataC,0);
resultT=analyze(dataT,0);

% for i=1:5
%     ScatterBox(resultC([i:5:length(resultC)],1),resultT([i:5:length(resultT)],1),{'Con','Tilt'},{'Stability (r)'},2)
% end
for i=1:5
    Group1(:,1,i)=resultC([i:5:length(resultC)],1);
    Group2(:,1,i)=resultT([i:5:length(resultT)],1);
end
c=1;
shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];
% CONTROL
y=[nanmean(Group1(:,c,1)),nanmean(Group1(:,c,2)),nanmean(Group1(:,c,3)),nanmean(Group1(:,c,4)),nanmean(Group1(:,c,5))];
SEM=[nanstd(Group1(:,c,1))/sqrt(size(Group1(:,c,1),1)),nanstd(Group1(:,c,2))/sqrt(size(Group1(:,c,2),1)),...
    nanstd(Group1(:,c,3))/sqrt(size(Group1(:,c,3),1)),nanstd(Group1(:,c,4))/sqrt(size(Group1(:,c,4),1)),...
    nanstd(Group1(:,c,5))/sqrt(size(Group1(:,c,5),1))];
h1=shadedErrorBar(1:5,y,SEM,'-k',0);
ylabel('Stability (r)');
xlabel('Sessions');
ax1=gca; set(ax1,'XTick',[1 2 3 4 5],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
hold on
% TILTED
y=[nanmean(Group2(:,c,1)),nanmean(Group2(:,c,2)),nanmean(Group2(:,c,3)),nanmean(Group2(:,c,4)),nanmean(Group2(:,c,5))];
SEM=[nanstd(Group2(:,c,1))/sqrt(size(Group2(:,c,1),1)),nanstd(Group2(:,c,2))/sqrt(size(Group2(:,c,2),1)),...
    nanstd(Group2(:,c,3))/sqrt(size(Group2(:,c,3),1)),nanstd(Group2(:,c,4))/sqrt(size(Group2(:,c,4),1)),...
    nanstd(Group2(:,c,5))/sqrt(size(Group2(:,c,5),1))];
h2=shadedErrorBar(1:5,y,SEM,'-r',1);
set(h2.mainLine,'Color',[1 0 0])


groups=[repmat({'Control'},length(Group1(:,1,1)),1);repmat({'Tilted'},length(Group2(:,1,1)),1)];
measureVal=[Group1(:,1,1),Group1(:,1,2),Group1(:,1,3),Group1(:,1,4),Group1(:,1,5);Group2(:,1,1),Group2(:,1,2),Group2(:,1,3),Group2(:,1,4),Group2(:,1,5)];
t=table(groups,measureVal(:,1),measureVal(:,2),measureVal(:,3),measureVal(:,4),measureVal(:,5),'VariableNames',{'species','meas1','meas2','meas3','meas4','meas5'});
Meas = table([1 2 3 4 5]','VariableNames',{'Measurements'});
% CREATE MODELS
% OVER ALL SESSIONS
rm1 = fitrm(t,'meas1-meas5~species','WithinDesign',Meas);
% SESSIONS 1 2 3 - rotation
rm2 = fitrm(t,'meas1-meas3~species','WithinDesign',Meas(1:3,1));
% SESSIONS 1 4 5 - dark
rm3 = fitrm(t,'meas1,meas4,meas5~species','WithinDesign',Meas([1,4,5],1));
% SESSIONS 1 3 5 - stability of standard sessions
rm4 = fitrm(t,'meas1,meas3,meas5~species','WithinDesign',Meas([1,3,5],1));

% tbl = mauchly(rm)
% ranovatbl = ranova(rm)
disp('Multivariate Tests *ALL SESSIONS*')
manovatbl=manova(rm1)
disp('Tests of Between-Subjects Effects *ALL SESSIONS*')
anovatbl=anova(rm1)

disp('Multivariate Tests *ROTATION*')
manovatbl=manova(rm2)
disp('Tests of Between-Subjects Effects *ROTATION*')
anova(rm2)

disp('Multivariate Tests *DARK*')
manovatbl=manova(rm3)
disp('Tests of Between-Subjects Effects *DARK*')
anova(rm3)

disp('Multivariate Tests *STABILITY*')
manovatbl=manova(rm4)
disp('Tests of Between-Subjects Effects *STABILITY*')
anova(rm4)
%          print(figure(1),'-bestfit','-dpdf', '-r600',[FigureLocation,filesep,'TiltedStabilitycell19s1_r0.53.pdf'])
v=1;
disp(['Control  S1: ',num2str(nanmean(Group1(:,v,1))),' ± ',num2str(nanstd(Group1(:,v,1)/sqrt(length(Group1(:,v,1))))),...
    ', S2: ',num2str(nanmean(Group1(:,v,2))),' ± ',num2str(nanstd(Group1(:,v,2)/sqrt(length(Group1(:,v,2))))),...
    ', S3: ',num2str(nanmean(Group1(:,v,3))),' ± ',num2str(nanstd(Group1(:,v,3)/sqrt(length(Group1(:,v,3))))),...
    ', S4: ',num2str(nanmean(Group1(:,v,4))),' ± ',num2str(nanstd(Group1(:,v,4)/sqrt(length(Group1(:,v,4))))),...
    ', S5: ',num2str(nanmean(Group1(:,v,5))),' ± ',num2str(nanstd(Group1(:,v,5)/sqrt(length(Group1(:,v,5)))))])

disp(['Tilted  S1: ',num2str(nanmean(Group2(:,v,1))),' ± ',num2str(nanstd(Group2(:,v,1)/sqrt(length(Group2(:,v,1))))),...
    ', S2: ',num2str(nanmean(Group2(:,v,2))),' ± ',num2str(nanstd(Group2(:,v,2)/sqrt(length(Group2(:,v,2))))),...
    ', S3: ',num2str(nanmean(Group2(:,v,3))),' ± ',num2str(nanstd(Group2(:,v,3)/sqrt(length(Group2(:,v,3))))),...
    ', S4: ',num2str(nanmean(Group2(:,v,4))),' ± ',num2str(nanstd(Group2(:,v,4)/sqrt(length(Group2(:,v,4))))),...
    ', S5: ',num2str(nanmean(Group2(:,v,5))),' ± ',num2str(nanstd(Group2(:,v,5)/sqrt(length(Group2(:,v,5)))))])


i=1
measureVal=[Group1(:,i,3),Group1(:,i,4)];
t=table(measureVal(:,1),measureVal(:,2),'VariableNames',{'meas1','meas2'});
Meas = table([1 2]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas2~1','WithinDesign',Meas);
disp('Control session 3 4')
ranova(rm)

measureVal=[Group2(:,i,3),Group2(:,i,4)];
t=table(measureVal(:,1),measureVal(:,2),'VariableNames',{'meas1','meas2'});
Meas = table([1 2]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas2~1','WithinDesign',Meas);
disp('Tilted session 3 4')
ranova(rm)

function R=analyze(data,plots)
cxsID=fieldnames(data);
for i=1:length(cxsID)
    Session=data.(cxsID{i});
    disp(cxsID{i})
    split=round(length(Session)/2);
    S1=Session(1:split,:);
    S2=Session(split+1:end,:);
    MinY=min(Session(:,3));MaxY=max(Session(:,3));MinX=min(Session(:,2));MaxX=max(Session(:,2));
    [map1]=binzzz(S1,S1(S1(:,5)==1,:),61,MaxX,MinX,MaxY,MinY);
    [map2]=binzzz(S2,S2(S2(:,5)==1,:),61,MaxX,MinX,MaxY,MinY);
    map1(isnan(map1))=0;
    map2(isnan(map2))=0;
    r = corrcoef(map1,map2);
    R(i,1)=r(2,1);
    
    %     if R(i,1)>.5 || R(i,1)<.25 && plots==1
    if plots==1
        subplot(1,2,1);set(gca,'color',[1 1 1]);
        upsamRateMap=PerfectCircRateMap(map1,0);
        imAlpha=ones(size(upsamRateMap));
        imAlpha(isnan(upsamRateMap))=0;
        imagesc(upsamRateMap,'AlphaData',imAlpha);
        box off;axis image;axis off;colormap jet;
        
        subplot(1,2,2);
        upsamRateMap=PerfectCircRateMap(map2,0);
        imAlpha=ones(size(upsamRateMap));
        imAlpha(isnan(upsamRateMap))=0;
        imagesc(upsamRateMap,'AlphaData',imAlpha);
        box off;axis image;axis off;colormap jet;
        disp(['r=',num2str(R(i,1))])
        test=1;
    end
end
end

function SmoothRateMap=binzzz(occMatrix,spks_VEL,track_length,MaxX,MinX,MaxY,MinY)
nBinsx = round(track_length/2.44); nBinsy = round(track_length/2.44);
edges{1} = linspace(MinY, MaxY, nBinsy+1);
edges{2} = linspace(MinX, MaxX, nBinsx+1);
Omatrix = hist3([occMatrix(:,3) occMatrix(:,2)],'Edges',edges);
Omatrix(end,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/60;
occ(occ<0.100)=0;
Smatrix = hist3([spks_VEL(:,3), spks_VEL(:,2)],'Edges',edges);
Smatrix(end,:) = [];
Smatrix(:,end) = [];
FilledRateMatrix = Smatrix./occ;
FilledRateMatrix(isinf(FilledRateMatrix))=0;
filtWidth = [5 5]; filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout');
end