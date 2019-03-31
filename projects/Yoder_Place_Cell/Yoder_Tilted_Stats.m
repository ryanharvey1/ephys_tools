% Yoder_Tilted_Stats
% clear ; 
close all ; clc ; clear 

%  CONTROL           TILTED
% [SR019,SR20,RH40] [SR011,SR012,RH48,RH26,SR21,SR24,SK107]

% NUMBER OF CELLS PER ANIMAL
% 
% CONTROL [# of place cells:54]
%   SR19:   191
%   SR20:   151
%   RH40:   80
%   total:  422
%   percent: 12.79%
% 
% TILTED [# of place cells:100] 
%   SR011:  290
%   SR012:  233
%   RH48:   212
%   RH26:   257
%   SR21:   104
%   SR24:   41
%   SK107:  67
%   total:  1204
%   percent: 8.3%
z=(.1279-.83)/sqrt(((54+100)/(422+1204))*(1-((54+100)/((422+1204)))*(1/422)+(1/1204)));

% [h,p, chi2stat,df] = prop_test([54 100] , [422 1204], 0);

addpath(('/Users/RyanHarvey/GoogleDrive/MatlabDir/RC_notBoxPlot'))
addpath(('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis'))
addpath(('/Users/RyanHarvey/GoogleDrive/MatlabDir/CircStat2012a'))
addpath(('/Users/RyanHarvey/GoogleDrive/MatlabDir/FMAToolbox'))
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/Yoder_Place_Cell')

FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';

load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/GW_Field_StatsPlaceCells_Tilted_Mice.mat')
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewData.mat')
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewData2.mat');
newdata2=ans;

allcells=newdata.textdata(2:end,1);
placecells=data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,1);
for i=1:length(placecells)
    newmeasures(i,:)=newdata.data(strcmpi(placecells(i),allcells),13:16);
end

allcells2=newdata2.textdata(2:end,1);
for i=1:length(placecells)
    newmeasures2(i,:)=newdata2.data(strcmpi(placecells(i),allcells),17:18);
end

data.data.Field_StatsPlaceCells0x2DTilted=[data.data.Field_StatsPlaceCells0x2DTilted,newmeasures,newmeasures2];

control=data.data.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,2),'Control'),:);
tilted=data.data.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(2:end,2),'Tilted'),:);

controlid=data.textdata.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(:,2),'Control'),1);
tiltedid=data.textdata.Field_StatsPlaceCells0x2DTilted(ismember(data.textdata.Field_StatsPlaceCells0x2DTilted(:,2),'Tilted'),1);

for i=1:5
    Group1id(:,:,i)=controlid(i:5:end,:);
    Group2id(:,:,i)=tiltedid(i:5:end,:);
end


for i=1:5
    Group1(:,:,i)=control(i:5:end,:);
    Group2(:,:,i)=tilted(i:5:end,:);
end

% ADD ECCENTRCITY
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/New_eccentricity.mat')
[r,c,d]=size(Group1);
for i=1:5
    Group1(:,c+1,i)=eccentricity.ec(:,1,i);
end
[r,c,d]=size(Group2);
for i=1:5
    Group2(:,c+1,i)=eccentricity.et(:,1,i);
end

% ADD Velocity
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/VelocityOverSessions.mat')
[r,c,d]=size(Group1);
for i=1:5
    Group1(:,c+1,i)=VelocityOverSessions.controlresult(:,1,i);
end
[r,c,d]=size(Group2);
for i=1:5
    Group2(:,c+1,i)=VelocityOverSessions.tiltedresult(:,1,i);
end


GroupNames={'Control','Tilted'};
VarNames=[data.textdata.Field_StatsPlaceCells0x2DTilted(1,3:end),newdata.textdata(1,end-3:end),newdata2.textdata(1,end-1:end)];
VarNames=regexprep(VarNames,'_','','emptymatch');

VarNames=[VarNames,'Eccentricity','Velocity'];


controldisp=Group1(:,15:16,2); Group1(:,15:16,:)=[];
tiltdisp=Group2(:,15:16,2); Group2(:,15:16,:)=[];
DispVarNames=VarNames(15:16); VarNames(15:16)=[];

%% INFIELD VS OUTFIELD RATIO
for i=1:5
    Group1(:,17,i)=Group1(:,15,i)./Group1(:,16,i);
    Group2(:,17,i)=Group2(:,15,i)./Group2(:,16,i);
end
% for i=1:5
%     Group1(:,17,i)=(Group1(:,15,i)-Group1(:,16,i))/(sum([Group1(:,15,i);Group1(:,16,i)]));
%     Group2(:,17,i)=(Group2(:,15,i)-Group2(:,16,i))/(sum([Group2(:,15,i);Group2(:,16,i)]));
% end
Group1(isinf(Group1(:,17)),17)=0;
Group2(isinf(Group2(:,17)),17)=0;
VarNames=[VarNames,'InfieldOutfieldRatio'];
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DISPLACEMENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% CDF WITH SHUFFLED DISTRIBUTION
load('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper/Displacement_Shuffled_correlations.mat')
cdfFig=figure;
cdfFig.Color=[1 1 1];
% PLOT SHUFFLE
[f,x] = ecdf(rSHUFF);
p4=plot(x,f);
set(p4,'LineWidth',5,'Color',[.5 .5 .5])
hold on
% PLOT CUT OFF
p4_1=plot(repmat(prctile(rSHUFF,95),2,1),[0;1],'LineStyle','--','Color',[.5 .5 .5],'LineWidth',5);
hold on
% PLOT CONTROL
[f,x] = ecdf(controldisp(:,2));
p4_2=plot(x,f);
set(p4_2,'LineWidth',5,'Color','k')
hold on
% PLOT TILTED
[f,x] = ecdf(tiltdisp(:,2));
p4_3=plot(x,f);
set(p4_3,'LineWidth',5,'Color','r')
box off
hold on
legend([p4 p4_1 p4_2 p4_3],{'Shuffled','95thPercentile','Control','Tilted'},'FontSize',20,'Location','best','box','off')
xlabel('Correlation');ylabel('Cumulative Frequency')
ax=gca;set(ax,'FontSize',20,'FontWeight','bold')
% print(cdfFig,'-bestfit', '-dpdf', '-r600',[FigureLocation,filesep,'cdf_Fig.pdf'])


% filter out low correlations
perBel95con=(sum(controldisp(:,2)<prctile(rSHUFF,95))/length(controldisp(:,2)))*100
perBel95tilt=(sum(tiltdisp(:,2)<prctile(rSHUFF,95))/length(tiltdisp(:,2)))*100
[h,p, chi2stat,df] = prop_test([sum(controldisp(:,2)<prctile(rSHUFF,95)) sum(tiltdisp(:,2)<prctile(rSHUFF,95))] , [length(controldisp(:,2)) length(tiltdisp(:,2))], 0)



controldisp(controldisp(:,2)<prctile(rSHUFF,95),:)=[];
tiltdisp(tiltdisp(:,2)<prctile(rSHUFF,95),:)=[];


% PIE CHART FOR AMOUNT OF CELLS THAT ROTATE
% don't rotate
norotcon=(sum(controldisp(:,1)>360-45 | controldisp(:,1)<45)/length(controldisp))*100
norottilt=(sum(tiltdisp(:,1)>360-45 | tiltdisp(:,1)<45)/length(tiltdisp))*100
% rotate with cue
rotcon=(sum(controldisp(:,1)>90-45 & controldisp(:,1)<90+45)/length(controldisp))*100
rottilt=(sum(tiltdisp(:,1)>90-45 & tiltdisp(:,1)<90+45)/length(tiltdisp))*100
% rotates not with cue 
otherrotcon=(sum(controldisp(:,1)>135 & controldisp(:,1)<315)/length(controldisp))*100
otherrottilt=(sum(tiltdisp(:,1)>135 & tiltdisp(:,1)<315)/length(tiltdisp))*100


% CHI SQUARE
x1=[ones(sum(controldisp(:,1)>360-45|controldisp(:,1)<45),1);ones(sum(controldisp(:,1)>90-45&controldisp(:,1)<90+45),1)+1;ones(sum(controldisp(:,1)>135 & controldisp(:,1)<315),1)+2];
x11=[ones(sum(tiltdisp(:,1)>360-45|tiltdisp(:,1)<45),1);ones(sum(tiltdisp(:,1)>90-45&tiltdisp(:,1)<90+45),1)+1;ones(sum(tiltdisp(:,1)>135 & tiltdisp(:,1)<315),1)+2];
x1_1=[x1;x11];
x2=[ones(length(x1),1);2*ones(length(x11),1)];
[table,chi2,p] = crosstab(x1_1,x2)

% test of porportion [fig 3 c pie charts]
% [h,p, chi2stat,df] = prop_test([88 61] , [100 100], 0)




pie_displacement_fig=figure; pie_displacement_fig.Color=[1 1 1];
ax=subplot(1,2,1);p1=pie([norotcon,rotcon,otherrotcon]);
title(ax,'Control');ax.FontSize=20;
hold on
ax=subplot(1,2,2);p2=pie([norottilt,rottilt,otherrottilt]);
title(ax,'Tilted');ax.FontSize=20;
print(pie_displacement_fig, '-dpng', '-r600',[FigureLocation,filesep,'pie_displacement_fig.png'])




controlstats=circ_stats(deg2rad(controldisp(:,1)))
tiltstats=circ_stats(deg2rad(tiltdisp(:,1)))


Polar_Fig=figure(4);
Polar_Fig.Color=[1 1 1];
p1=polarhistogram(deg2rad(controldisp(:,1)),60,'FaceColor','k','FaceAlpha',0.8,'EdgeColor','w','Normalization','probability'); hold on
p2=polarhistogram(deg2rad(tiltdisp(:,1)),60,'FaceColor','red','FaceAlpha',.5,'EdgeColor','w','Normalization','probability');
ax=gca; ax.ThetaTick=[0,90,180,270]; ax.FontWeight='bold'; ax.FontSize=20; ax.RAxisLocation=45;ax.GridAlpha=.5;ax.GridColor='k';
ax.RTick=[round(linspace(0,.16,4),2)];
% print(Polar_Fig,'-bestfit', '-dpdf', '-r300',[FigureLocation,filesep,'Polar_Fig.pdf'])

Polar_Fig1=figure(4);
Polar_Fig1.Color=[1 1 1];
p1=polarhistogram(deg2rad(controldisp(:,1)),60,'FaceColor','k','FaceAlpha',0.8,'EdgeColor','w','Normalization','probability'); hold on
ax=gca; ax.ThetaTick=[0,90,180,270]; ax.FontWeight='bold'; ax.FontSize=20; ax.RAxisLocation=45;ax.GridAlpha=.5;ax.GridColor='k'; ax.RLim=[0 0.16]
ax.RTick=[round(linspace(0,.16,4),2)];

Polar_Fig2=figure(5);
Polar_Fig2.Color=[1 1 1];
p2=polarhistogram(deg2rad(tiltdisp(:,1)),60,'FaceColor','red','FaceAlpha',.8,'EdgeColor','w','Normalization','probability');
ax=gca; ax.ThetaTick=[0,90,180,270]; ax.FontWeight='bold'; ax.FontSize=20; ax.RAxisLocation=45;ax.GridAlpha=.5;ax.GridColor='k'; ax.RLim=[0 0.16]
ax.RTick=[round(linspace(0,.16,4),2)];

print(Polar_Fig1,'-bestfit', '-dpdf', '-r600',[FigureLocation,filesep,'Polar_Fig_control.pdf'])
print(Polar_Fig2,'-bestfit', '-dpdf', '-r600',[FigureLocation,filesep,'Polar_Fig_tilted.pdf'])



% Wilcoxon rank sum test
%  SUBTRACK EACH ANGLE FROM MEAN AND THEN RUN TEST
[mu, ~, ~] = circ_mean(deg2rad(controldisp(:,1))); x=controldisp(:,1)-rad2deg(mu);
[mu, ~, ~] = circ_mean(deg2rad(tiltdisp(:,1))); y=tiltdisp(:,1)-rad2deg(mu);
[p,h,stats] = ranksum(x,y)


% resultant vector
rcontrol = circ_r(deg2rad(controldisp(:,1)))
rtilt = circ_r(deg2rad(tiltdisp(:,1)))
Rtogether=circ_r([[deg2rad(controldisp(:,1));[deg2rad(tiltdisp(:,1))]]])

[pval, f] = circ_ktest(deg2rad(controldisp(:,1)),deg2rad(tiltdisp(:,1)))
[pval, table] =circ_wwtest(deg2rad(controldisp(:,1)),deg2rad(tiltdisp(:,1)))


[p,U2_obs,U2_H0]=watsons_U2_perm_test(deg2rad(controldisp(:,1)),deg2rad(tiltdisp(:,1)),10)


[p alpha] = circ_vmpdf(deg2rad(controldisp(:,1)), 90, 1)


group=[ones(length(controldisp(:,1)),1); ones(length(tiltdisp(:,1)),1)+1];
% [p,F] = CircularANOVA(deg2rad(controldisp(:,1)),deg2rad(tiltdisp(:,1)),group,'lr')

[h,p] = ConcentrationTest([deg2rad(controldisp(:,1));deg2rad(tiltdisp(:,1))],group)

% Watson's Goodness-of-fit tests
[U2,p] = watson1962(controldisp(:,1),tiltdisp(:,1))
% Kruskal-Wallis
[pval,med,P] = circ_cmtest(deg2rad(controldisp(:,1)),deg2rad(tiltdisp(:,1)))

PhaseStats=[strcat('Displacement',': ','Kruskal-Wallis Statistic:',num2str(P),'  Shared Medium:',num2str(rad2deg(med)),'  PVal:',num2str(pval)),...
    strcat('Watson''s Goodness of fit: ','  U2:',num2str(U2),'  Pvalue:',num2str(p))];

for c=[1,2,3,7,8,9,10,12,13,14]
    shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 255 421];
    % CONTROL
    y=[nanmean(Group1(:,c,1)),nanmean(Group1(:,c,2))];
    SEM=[nanstd(Group1(:,c,1))/sqrt(size(Group1(:,c,1),1)),nanstd(Group1(:,c,2))/sqrt(size(Group1(:,c,2),1))];
    h1=shadedErrorBar(1:2,y,SEM,'-k',0);
    ylabel(VarNames(c));
    xlabel('Session')
    ax1=gca; set(ax1,'XTick',[1 2],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
    
    hold on
    
    % TILTED
    y=[nanmean(Group2(:,c,1)),nanmean(Group2(:,c,2))];
    SEM=[nanstd(Group2(:,c,1))/sqrt(size(Group2(:,c,1),1)),nanstd(Group2(:,c,2))/sqrt(size(Group2(:,c,2),1))];
    h2=shadedErrorBar(1:2,y,SEM,'-r',1);
    set(h2.mainLine,'Color',[1 0 0])
    
%     print(shadded_line_fig, '-dpdf', '-r300',{strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')})
%     print(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig.pdf'),'-dpdf','-r300','-bestfit')
    print(shadded_line_fig,'-bestfit', '-dpdf', '-r600',char(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig_3')))

    close all
end



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[ AllStats ] = ScatterBox(Group1(:,:,1),Group2(:,:,1),GroupNames,VarNames,2)
print(figure(1),'-bestfit', '-dpdf', '-r300',[FigureLocation,filesep,'Field2wall_Scatter_Dot.pdf'])

pause
[ AllStats ] = ScatterBox(controldisp(:,2),tiltdisp(:,2),GroupNames,DispVarNames(2),2)


% STAIR HIST
% Controldisp=histcounts(Group1(:,10,1),22);
% Tilteddisp=histcounts(Group2(:,10,1),22);
%   Controldisp= (Controldisp - min(Controldisp)) / ( max(Controldisp) - min(Controldisp) );
%   Tilteddisp= (Tilteddisp - min(Tilteddisp)) / ( max(Tilteddisp) - min(Tilteddisp) );
% h1=bar(Controldisp); hold on; h2=bar(Tilteddisp);

% ECDF
data1=Group1(:,13,1);
data2=Group2(:,13,1);
[p,h,stats] = ranksum(data1,data2)

[f,x] = ecdf(data1);
p4=plot(x,f); x1=x;
set(p4,'LineWidth',10,'Color','k')
hold on
[f,x] = ecdf(data2);
p4_2=plot(x,f); x2=x;
set(p4_2,'LineWidth',10,'Color','r')
box off
legend([p4 p4_2],{'Control','PAE'},'FontSize',20,'Location','best')
xlabel(VarNames(13));ylabel('Cumulative Frequency')

% PIE CHART BORDER SCORE
pie_border_fig=figure; pie_border_fig.Color=[1 1 1];
ax=subplot(1,2,1);p1=pie([sum(data1>0)/length(data1),sum(data1<0)/length(data1)]);
title(ax,'Control');ax.FontSize=20;
hold on
ax=subplot(1,2,2);p2=pie([sum(data2>0)/length(data2),sum(data2<0)/length(data2)]);
title(ax,'Tilted');ax.FontSize=20;
print(pie_border_fig,'-bestfit', '-dpdf', '-r300',[FigureLocation,filesep,'pie_border_fig.pdf'])



% REMOVE UNDERSCORES FROM VAR NAMES
VarNames=regexprep(VarNames,'_','','emptymatch');
%%
% SHADED LINE GRAPH
% shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];

% cc=1;
for c=[1,2,3,7,8,9,10,12,13,14,15,16,17]
    shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];
    % CONTROL
    y=[nanmean(Group1(:,c,1)),nanmean(Group1(:,c,2)),nanmean(Group1(:,c,3)),nanmean(Group1(:,c,4)),nanmean(Group1(:,c,5))];
    SEM=[nanstd(Group1(:,c,1))/sqrt(size(Group1(:,c,1),1)),nanstd(Group1(:,c,2))/sqrt(size(Group1(:,c,2),1)),...
        nanstd(Group1(:,c,3))/sqrt(size(Group1(:,c,3),1)),nanstd(Group1(:,c,4))/sqrt(size(Group1(:,c,4),1)),...
        nanstd(Group1(:,c,5))/sqrt(size(Group1(:,c,5),1))];
    h1=shadedErrorBar(1:5,y,SEM,'-k',0);
    ylabel(VarNames(c));
    ax1=gca; set(ax1,'XTick',[1 2 3 4 5],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
    
    hold on
    
    % TILTED
    y=[nanmean(Group2(:,c,1)),nanmean(Group2(:,c,2)),nanmean(Group2(:,c,3)),nanmean(Group2(:,c,4)),nanmean(Group2(:,c,5))];
    SEM=[nanstd(Group2(:,c,1))/sqrt(size(Group2(:,c,1),1)),nanstd(Group2(:,c,2))/sqrt(size(Group2(:,c,2),1)),...
        nanstd(Group2(:,c,3))/sqrt(size(Group2(:,c,3),1)),nanstd(Group2(:,c,4))/sqrt(size(Group2(:,c,4),1)),...
        nanstd(Group2(:,c,5))/sqrt(size(Group2(:,c,5),1))];
    h2=shadedErrorBar(1:5,y,SEM,'-r',1);
    set(h2.mainLine,'Color',[1 0 0])
    
%     print(shadded_line_fig, '-dpdf', '-r300',{strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')})
%     print(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig.pdf'),'-dpdf','-r300','-bestfit')
    print(shadded_line_fig,'-bestfit', '-dpdf', '-r600',char(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')))

    close all
end

%%
% SHADDED LINE GRAPHS FOR DARK SESSION
for c=[1,2,3,7,8,9,10,12,13,14,15,16]
    shadded_line_dark_fig=figure; shadded_line_dark_fig.Color=[1 1 1];shadded_line_dark_fig.OuterPosition=[977 418 428 410];
    % CONTROL
    y=[nanmean(Group1(:,c,3)),nanmean(Group1(:,c,4)),nanmean(Group1(:,c,5))];
    SEM=[nanstd(Group1(:,c,3))/sqrt(size(Group1(:,c,3),1)),nanstd(Group1(:,c,4))/sqrt(size(Group1(:,c,4),1)),nanstd(Group1(:,c,5))/sqrt(size(Group1(:,c,5),1))];
    
    h1=shadedErrorBar(3:5,y,SEM,'-k',0);
    ylabel(VarNames(c));
    xlabel('Session')
    ax1=gca; set(ax1,'XTick',[3 4 5],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
    
    hold on
    
    % TILTED
    y=[nanmean(Group2(:,c,3)),nanmean(Group2(:,c,4)),nanmean(Group2(:,c,5))];
    SEM=[nanstd(Group2(:,c,3))/sqrt(size(Group2(:,c,3),1)),nanstd(Group2(:,c,4))/sqrt(size(Group2(:,c,4),1)),nanstd(Group2(:,c,5))/sqrt(size(Group2(:,c,5),1))];
    h2=shadedErrorBar(3:5,y,SEM,'-r',1);
    set(h2.mainLine,'Color',[1 0 0])
    
%     print(shadded_line_fig, '-dpdf', '-r300',{strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')})
%     print(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig.pdf'),'-dpdf','-r300','-bestfit')
    print(shadded_line_dark_fig, '-dpdf', '-r600',char(strcat(FigureLocation,filesep,VarNames(c),'shadded_line_dark_fig')))

    close all
end

% SHADDED LINE GRAPHS FOR stability (sess 1 3 5)
for c=[1,2,3,7,8,9,10,12,13,14]
    shadded_line_stability_fig=figure; shadded_line_stability_fig.Color=[1 1 1];shadded_line_stability_fig.OuterPosition=[977 418 428 410];
    % CONTROL
    y=[nanmean(Group1(:,c,1)),nanmean(Group1(:,c,3)),nanmean(Group1(:,c,5))];
    SEM=[nanstd(Group1(:,c,1))/sqrt(size(Group1(:,c,1),1)),nanstd(Group1(:,c,3))/sqrt(size(Group1(:,c,3),1)),nanstd(Group1(:,c,5))/sqrt(size(Group1(:,c,5),1))];
    h1=shadedErrorBar([1,3,5],y,SEM,'-k',0);
    ylabel(VarNames(c));
    xlabel('Session')
    ax1=gca; set(ax1,'XTick',[1 3 5],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
    
    hold on
    
    % TILTED
    y=[nanmean(Group2(:,c,1)),nanmean(Group2(:,c,3)),nanmean(Group2(:,c,5))];
    SEM=[nanstd(Group2(:,c,1))/sqrt(size(Group2(:,c,1),1)),nanstd(Group2(:,c,3))/sqrt(size(Group2(:,c,3),1)),nanstd(Group2(:,c,5))/sqrt(size(Group2(:,c,5),1))];
    h2=shadedErrorBar([1,3,5],y,SEM,'-r',1);
    set(h2.mainLine,'Color',[1 0 0])
    
    print(shadded_line_stability_fig, '-dpdf', '-r600',char(strcat(FigureLocation,filesep,VarNames(c),'shadded_line_stability_fig')))

    close all
end
% SIMPLE EFFECTS FOR STABILITY
% session 1 to 3
i=9;
groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
% groups=[ones(length(Group1(:,i,1)),1);ones(length(Group2(:,i,1)),1)+1];
measureVal=[Group1(:,i,3),Group1(:,i,4);Group2(:,i,3),Group2(:,i,4)];
t=table(groups,measureVal(:,1),measureVal(:,2),'VariableNames',{'species','meas1','meas2'});
Meas = table([1 2]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas2~species','WithinDesign',Meas);
disp('session 3 to 4')
disp(['Control S3 mean: ', num2str(nanmean(Group1(:,i,3))), '   Control S4 mean: ', num2str(nanmean(Group1(:,i,4)))])
disp(['Tilted S3 mean:  ', num2str(nanmean(Group2(:,i,3))), '   Tilted S4 mean:  ', num2str(nanmean(Group2(:,i,4)))])

manova(rm)

% session 1 to 5
groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
measureVal=[Group1(:,i,3),Group1(:,i,5);Group2(:,i,3),Group2(:,i,5)];
t=table(groups,measureVal(:,1),measureVal(:,2),'VariableNames',{'species','meas1','meas2'});
Meas = table([1 2]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas2~species','WithinDesign',Meas);
disp('session 3 to 5')
disp(['Control S3 mean: ', num2str(nanmean(Group1(:,i,3))), '   Control S5 mean: ', num2str(nanmean(Group1(:,i,5)))])
disp(['Tilted S3 mean:  ', num2str(nanmean(Group2(:,i,3))), '   Tilted S5 mean:  ', num2str(nanmean(Group2(:,i,5)))])

manova(rm)

% session 3 to 5
groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
measureVal=[Group1(:,i,4),Group1(:,i,5);Group2(:,i,4),Group2(:,i,5)];
t=table(groups,measureVal(:,1),measureVal(:,2),'VariableNames',{'species','meas1','meas2'});
Meas = table([1 2]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas2~species','WithinDesign',Meas);
disp('session 4 to 5')
disp(['Control S4 mean: ', num2str(nanmean(Group1(:,i,4))), '   Control S5 mean: ', num2str(nanmean(Group1(:,i,5)))])
disp(['Tilted S4 mean:  ', num2str(nanmean(Group2(:,i,4))), '   Tilted S5 mean:  ', num2str(nanmean(Group2(:,i,5)))])

manova(rm)





%%
clc
% session 1 to 5
for i=1:length(VarNames)
groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
measureVal=[Group1(:,i,1),Group1(:,i,3),Group1(:,i,5)];
t=table(measureVal(:,1),measureVal(:,2),measureVal(:,3),'VariableNames',{'meas1','meas2','meas3'});
Meas = table([1 2 3]','VariableNames',{'Measurements'});
rm = fitrm(t,'meas1-meas3~1','WithinDesign',Meas);
disp(strcat('--------------------------------------------------',VarNames(i),'--------------------------------------------------'))
disp('session 1 3 5')
ranova(rm)
end

% % REPEATED MEASURES MODEL
% for ii=1:length(VarNames)
%     results(ii).varnames = VarNames(ii);
% end
%%
for i=1:length(VarNames)
groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
measureVal=[Group1(:,i,1),Group1(:,i,2),Group1(:,i,3),Group1(:,i,4),Group1(:,i,5);Group2(:,i,1),Group2(:,i,2),Group2(:,i,3),Group2(:,i,4),Group2(:,i,5)];
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
disp(strcat('--------------------------------------------------',VarNames(i),'--------------------------------------------------'))
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
end
%%
% SIMPLE EFFECTS FOR ABOVE TESTS
i=9;
groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
measureVal=[Group1(:,i,1),Group1(:,i,2);Group2(:,i,1),Group2(:,i,2)];
t=table(groups,measureVal(:,1),measureVal(:,2),'VariableNames',{'species','meas1','meas2'});
Meas = table([1 2]','VariableNames',{'Measurements'});

rm = fitrm(t,'meas1-meas2~species','WithinDesign',Meas);

anova(rm)

i=9;
groups=data.textdata.Session1; groups(1,:)=[]; groups=groups(~cellfun(@isempty,groups));
measureVal=[Group1(:,i,1),Group1(:,i,3);Group2(:,i,1),Group2(:,i,3)];
t=table(groups,measureVal(:,1),measureVal(:,2),'VariableNames',{'species','meas1','meas2'});
Meas = table([1 2]','VariableNames',{'Measurements'});

rm = fitrm(t,'meas1-meas2~species','WithinDesign',Meas);

anova(rm)




% COMPILE FOR SPSS OUTPUT

vars=repmat(VarNames,[1 5]); 
for iii=1:14
    ii=1;
    for i=iii:14:70
        vars{i}=strcat(vars{i},num2str(ii));ii=ii+1;
    end
end


['Group',repmat(VarNames,[1 5])]
SPSS=[[ones(size(Group1(:,:,1),1),1);ones(size(Group2(:,:,1),1),1)+1],[[Group1(:,:,1),Group1(:,:,2),Group1(:,:,3),Group1(:,:,4),Group1(:,:,5)];[Group2(:,:,1),Group2(:,:,2),Group2(:,:,3),Group2(:,:,4),Group2(:,:,5)]]];
table(SPSS,'VariableNames',['Group',vars]);



% BETWEEN SUBJECTS COLAPSE FOR FIGURE 1

controlDat=[Group1(:,:,1);Group1(:,:,2);Group1(:,:,3);Group1(:,:,4);Group1(:,:,5)];
TiltedDat=[Group2(:,:,1);Group2(:,:,2);Group2(:,:,3);Group2(:,:,4);Group2(:,:,5)];

% remove error outlier in coherence
controlDat(controlDat(:,9)>1,9)=controlDat((controlDat(:,9)>1),9)*0.1;


% [ AllStats ] = ScatterBox(controlDat(:,[1,9,14]),TiltedDat(:,[1,9,14]),GroupNames,VarNames(:,[1,9,14]),1)
% print(figure(1),'-bestfit', '-dpdf', '-r600',[FigureLocation,filesep,'between_Subs_Fig1.pdf'])

[ AllStats ] = ScatterBox(controlDat(:,1),TiltedDat(:,1),GroupNames,VarNames(:,1),1)
print(figure(1),'-bestfit', '-dpdf', '-r600',[FigureLocation,filesep,'between_Subs_PeakFR_Fig1.pdf'])

[ AllStats ] = ScatterBox(controlDat(:,9),TiltedDat(:,9),GroupNames,VarNames(:,9),1)
print(figure(2),'-bestfit', '-dpdf', '-r600',[FigureLocation,filesep,'between_Subs_Coherence_Fig1.pdf'])

[ AllStats ] = ScatterBox(controlDat(:,14),TiltedDat(:,14),GroupNames,VarNames(:,14),1)
print(figure(3),'-bestfit', '-dpdf', '-r600',[FigureLocation,filesep,'between_Subs_FieldWidth_Fig1.pdf'])





data1=controlDat(:,14);
data2=TiltedDat(:,14);

[f,x] = ecdf(data1);
p4=plot(x,f); x1=x;
set(p4,'LineWidth',10,'Color','k')
hold on
[f,x] = ecdf(data2);
p4_2=plot(x,f); x2=x;
set(p4_2,'LineWidth',10,'Color','r')
box off
legend([p4 p4_2],{'Control','Tilted'},'FontSize',20,'Location','best')
xlabel(VarNames(1));ylabel('Cumulative Frequency')




%% 
% LOAD SS3D TT FILE
% fn='/Users/RyanHarvey/Desktop/TT4 - MARGINAL2_2.ntt';
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analyses/spikeCode/MEX'));
fn= '/Volumes/Ryan_4TB/Place_Cell_Data/Place Cells - Tilted Mice/2013-10-28_13-51-26 - RH026 - 1 Place cell, analyzed/TT2- RH.ntt';
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] =Nlx2MatSpike_v3(fn,[1 1 1 1 1],1,1,[]);
cell1=Samples(:,:,CellNumbers==1);
cell2=Samples(:,:,CellNumbers==2);
cell3=Samples(:,:,CellNumbers==3);
cell4=Samples(:,:,CellNumbers==4);
cell5=Samples(:,:,CellNumbers==5);
cell6=Samples(:,:,CellNumbers==6);
cell7=Samples(:,:,CellNumbers==7);
cell8=Samples(:,:,CellNumbers==8);

[r,c,d]=size(cell7);
waveforms=figure; waveforms.Color=[1 1 1];
for ii=1:4
    subplot(2,2,ii)
    for i=1:round(d/4)
        plot(cell7(:,ii,i))
        hold on
    end
end

print(figure(1), '-dpdf', '-r600',[FigureLocation,filesep,'WaveformFinal.pdf'])

%%

% LOCOMOTION
load('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis/Yoder_Place_Cell/Locomotion','Controlresults','Tiltedresults');
[ AllStats ] = ScatterBox(Controlresults,Tiltedresults,{'Control','Tilted'},{'Head Angle per Sec','Distance','Mean Velocity','STDVelocity','medVelocity'},2)

% correlations
figure;
subplot(2,2,1)
scatter(Controlresults(:,3),mean(Group1(:,1,:),3))
xlabel('speed');ylabel('peak rate')
title(['Control Correlation: ',num2str(corr(mean(Group1(:,1,:),3),Controlresults(:,3)))])

subplot(2,2,2)
scatter(Controlresults(:,3),mean(Group1(:,14,:),3))
xlabel('speed');ylabel('field width')
title(['Control Correlation: ',num2str(corr(mean(Group1(:,14,:),3),Controlresults(:,3)))])

subplot(2,2,3)
scatter(Tiltedresults(:,3),mean(Group2(:,1,:),3))
xlabel('speed');ylabel('peak rate')
title(['Tilted Correlation: ',num2str(corr(mean(Group1(:,1,:),3),Controlresults(:,3)))])

subplot(2,2,4)
scatter(Tiltedresults(:,3),mean(Group2(:,14,:),3))
xlabel('speed');ylabel('field width')
title(['Tilted Correlation: ',num2str(corr(mean(Group1(:,14,:),3),Controlresults(:,3)))])

ControlMean=mean(Controlresults(:,3)) 
ControlSEM=std(Controlresults(:,3))/sqrt(length(Controlresults(:,3)))
TiltedMean=mean(Tiltedresults(:,3)) 
TiltedSEM=std(Tiltedresults(:,3))/sqrt(length(Tiltedresults(:,3)))
% 
% 
% FIELDWIDTH TO COHERENCE CORRELATIONS
figure;
subplot(1,2,1)
Data=[Group1(:,14,1),Group1(:,9,1)];
Data(logical(sum(isnan(Data),2)),:)=[];

scatter(Data(:,1),Data(:,2))
xlabel('Field Width');ylabel('Coherence')
title(['Control Correlation: ',num2str(corr(Data(:,1),Data(:,2)))])
lsline

clear Data

subplot(1,2,2)
Data=[Group2(:,14,1),Group2(:,9,1)];
Data(logical(sum(isnan(Data),2)),:)=[];

scatter(Data(:,1),Data(:,2))
xlabel('Field Width');ylabel('Coherence')
title(['Tilted Correlation: ',num2str(corr(Data(:,1),Data(:,2)))])
lsline


%% MEANS AND SEM
% COLAPSED
for v=1:length(VarNames)
disp(VarNames(v))
ControlMean=nanmean([Group1(:,v,1);Group1(:,v,2);Group1(:,v,3);Group1(:,v,4);Group1(:,v,5)])
ControlSEM=nanstd([Group1(:,v,1);Group1(:,v,2);Group1(:,v,3);Group1(:,v,4);Group1(:,v,5)])/sqrt(length([Group1(:,v,1);Group1(:,v,2);Group1(:,v,3);Group1(:,v,4);Group1(:,v,5)]))

TiltedMean=nanmean([Group2(:,v,1);Group2(:,v,2);Group2(:,v,3);Group2(:,v,4);Group2(:,v,5)])
TiltedSEM=nanstd([Group2(:,v,1);Group2(:,v,2);Group2(:,v,3);Group2(:,v,4);Group2(:,v,5)])/sqrt(length([Group2(:,v,1);Group2(:,v,2);Group2(:,v,3);Group2(:,v,4);Group2(:,v,5)]))
end
%% BY SESSION
clc
% compute stats for angular data from session 2 first
displacementMeanCon=rad2deg(circ_mean(deg2rad(Group1(:,15,2))))
displacementMeanTilt=rad2deg(circ_mean(deg2rad(Group2(:,15,2))))

displacementSEMCon=rad2deg(circ_std(deg2rad(Group1(:,15,2)/sqrt(length(Group1(:,15,2))))))
displacementSEMPAE=rad2deg(circ_std(deg2rad(Group2(:,15,2)/sqrt(length(Group2(:,15,2))))))

% convert bins to cm 
for ii=1:5
    for i=[5,6,10,11,12,14]
        Group1(:,i,ii)=Group1(:,i,ii)*2.44; 
        Group2(:,i,ii)=Group2(:,i,ii)*2.44; 
    end
end

for v=1:length(VarNames)
    disp(VarNames(v))
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
end


%% Finding examples of low active to field width cells
active=normc(Group2(:,3,1));
width=normc(Group2(:,14,1));

[m,i]=min(Group2(:,9,1))

%%
close all
for i=1:5
    figure(i);scatter(Group1(:,1,i),(Group1(:,10,2)));hold on; scatter(Group2(:,1,i),(Group2(:,10,i)));
    
    [rho,pval]=corr(Group1(:,1,i),Group1(:,10,2));
    disp(['Control Cor for session ',num2str(i),': ',num2str(rho),' ',num2str(pval)])
    [rho,pval]=corr(Group2(:,1,i),Group2(:,10,2));
    disp(['Tilted Cor for session ',num2str(i),': ',num2str(rho),' ',num2str(pval)])
end



%% FIELD TO WALL / ECCENTRICITY CORRELATIONS
clear corC corT
for i=1:5
    [r,p]=corr(Group1(:,10,i),Group1(:,17,i),'Rows','complete');
    corC(i,:)=[r,p];
    [r,p]=corr(Group2(:,10,i),Group2(:,17,i),'Rows','complete');
    corT(i,:)=[r,p];
end
corC
corT

% for i=1:5
% F2W=Group1(:,11,i);E=Group1(:,17,i);
% temp=[F2W,E]; temp(isnan(F2W))=[];temp(isnan(E))=[];
% F2W=temp(:,1);
% E=temp(:,2);
% 
%     pc(i,:) = polyfit(F2W,E,1);
%     pt(i,:) = polyfit(F2W,E,1);
% end
% pc
% pt


for i=1:5
    b1=Group1(:,10,i)\Group1(:,17,i);
    yCalc1 = b1*Group1(:,10,i);
    X = [ones(length(Group1(:,10,i)),1) Group1(:,10,i)];
    b = X\Group1(:,17,i);
    yCalc2 = X*b;
    Rsq1c(i) = 1 - nansum((Group1(:,17,i) - yCalc1).^2)/nansum((Group1(:,17,i) - nanmean(Group1(:,17,i))).^2);
    Rsq2c(i) = 1 - nansum((Group1(:,17,i) - yCalc2).^2)/nansum((Group1(:,17,i) - nanmean(Group1(:,17,i))).^2);

    b1=Group2(:,11,i)\Group2(:,17,i);
    yCalc1 = b1*Group2(:,11,i);
    X = [ones(length(Group2(:,11,i)),1) Group2(:,11,i)];
    b = X\Group2(:,17,i);
    yCalc2 = X*b;
    Rsq1t(i) = 1 - nansum((Group2(:,17,i) - yCalc1).^2)/nansum((Group2(:,17,i) - nanmean(Group2(:,17,i))).^2);
    Rsq2t(i) = 1 - nansum((Group2(:,17,i) - yCalc2).^2)/nansum((Group2(:,17,i) - nanmean(Group2(:,17,i))).^2);
end
Rsq1c
Rsq2c
Rsq1t
Rsq2t

i=1;
lm = fitlm(Group1(:,10,i),Group1(:,17,i),'linear')
%% Eccentricity
split=median([Group1(:,10,1);Group2(:,10,1)]);

F2W=Group1(:,10,1);
E=Group1(:,17,1);

temp=[F2W,E]; 
temp(isnan(temp(:,1)),:)=[];
temp(isnan(temp(:,2)),:)=[];

awayc=temp(temp(:,1)>split,:);
close2wallc=temp(temp(:,1)<split,:);

clear temp F2W E

F2W=Group2(:,10,1);
E=Group2(:,17,1);

temp=[F2W,E]; 
temp(isnan(temp(:,1)),:)=[];
temp(isnan(temp(:,2)),:)=[];

awayt=temp(temp(:,1)>split,:);
close2wallt=temp(temp(:,1)<split,:);



[ AllStats ] = ScatterBox(awayc(:,2),awayt(:,2),{'Control','Tilted'},{'Eccentricity Away'},2)
[ AllStats ] = ScatterBox(close2wallc(:,2),close2wallt(:,2),{'Control','Tilted'},{'Eccentricity Close'},2)

[ AllStats ] = ScatterBox(awayc(:,2),close2wallc(:,2),{'Control','Tilted'},{'Control Eccentricity'},2)
[ AllStats ] = ScatterBox(awayt(:,2),close2wallt(:,2),{'Control','Tilted'},{'Tilted Eccentricity'},2)

[ AllStats ] = ScatterBox(Group1(:,17,1),Group2(:,17,1),{'Control','Tilted'},{'Between Group Eccentricity'},2)





%% ____________________________________________________K-Means Clusting _______________________________________________________________
rng(14,'twister');

meas=[Group1(:,:,1);Group2(:,:,1)];
groups=[repmat({'Group1'},length(Group1(:,:,1)),1);repmat({'Group2'},length(Group2(:,:,1)),1)];

groups(sum(isnan(meas),2)>0,:)=[];

meas(sum(isnan(meas),2)>0,:)=[];

%%
[cidx,cmeans] = kmeans(meas,3,'dist','sqeuclidean','display','iter');
[silh2,h] = silhouette(meas,cidx,'sqeuclidean');
%%

ptsymb = {'bs','r^','md','go','c+'};
for i = 1:3
    clust = find(cidx==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
    hold on
end
plot3(cmeans(:,1),cmeans(:,2),cmeans(:,3),'ko');
plot3(cmeans(:,1),cmeans(:,2),cmeans(:,3),'kx');
hold off
xlabel('Sepal Length');
ylabel('Sepal Width');
zlabel('Petal Length');
view(-137,10);
grid on

%%
eucD = pdist(meas,'euclidean');
clustTreeEuc = linkage(eucD,'average');
cophenet(clustTreeEuc,eucD)

figure(1)
[h,nodes] = dendrogram(clustTreeEuc,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];

cosD = pdist(meas,'cosine');
clustTreeCos = linkage(cosD,'average');
cophenet(clustTreeCos,cosD)

figure(2)
[h,nodes] = dendrogram(clustTreeCos,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];

figure(3)
[h,nodes] = dendrogram(clustTreeCos,12);

%%
[sum(ismember(nodes,[4 11 2 3])) sum(ismember(nodes,[8 12 5])) ...
                  sum(ismember(nodes,[6 7 9 10])) sum(nodes==1)]

%%
figure(4)
hidx = cluster(clustTreeCos,'criterion','distance','cutoff',.006);
for i = 1:5
    clust = find(hidx==i);
    plot3(meas(clust,1),meas(clust,2),meas(clust,3),ptsymb{i});
    hold on
end
hold off
xlabel('Sepal Length');
ylabel('Sepal Width');
zlabel('Petal Length');
view(-137,10);
grid on

%%
figure(5)
clustTreeSng = linkage(eucD,'single');
[h,nodes] = dendrogram(clustTreeSng,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];


%%
c = cluster(clustTreeSng,'maxclust',2);
[table,chi2,p]=crosstab(c,groups)



%%____________________________________________________________________________________________________________________________________

% FIELDTRIP
% [spike] = ft_read_spike(fn);
% 
% 
% cfg              = [];
% cfg.spikechannel = {'sig002a_wf', 'sig003a_wf'}; % select only the two single units
% spike = ft_spike_select(cfg, spike);
% 
% 
% cfg             = [];
% cfg.fsample     = 40000;
% cfg.interpolate = 1; % keep the density of samples as is
% [wave, spikeCleaned] = ft_spike_waveform(cfg,spike);
% 
% 
% for k = [1 2]
%   figure, 
%   sl = squeeze(wave.dof(k,:,:))>1000; % only keep samples with enough spikes
%   plot(wave.time(sl), squeeze(wave.avg(k,:,sl)),'k') % factor 10^6 to get microseconds
%   hold on
%  
%   % plot the standard deviation
%   plot(wave.time(sl), squeeze(wave.avg(k,:,sl))+sqrt(squeeze(wave.var(k,:,sl))),'k--') 
%   plot(wave.time(sl), squeeze(wave.avg(k,:,sl))-sqrt(squeeze(wave.var(k,:,sl))),'k--')
%  
%   axis tight
%   set(gca,'Box', 'off') 
%   xlabel('time')
%   ylabel('normalized voltage')
% end
%%
% vars=['Group',repmat(VarNames,[1 5])];

% fields=fieldnames(data.data);
% for i=2:6
%     measures=data.data.(fields{i});
%     split=sum(isnan(measures),2);
%     [~,I]=max(split);
%     data.data.(fields{i})=[data.data.(fields{i}),[newmeasures(1:I-1,:);NaN(1,size(newmeasures,2));newmeasures(I+1:size(data.data.(fields{i}),1),:)]];
%     newmeasures(1:size(data.data.(fields{i})-1,1),:)=[];
% end
% SORT OUT NEW DATA
% fields=fieldnames(data.data);
% ii=1;
% for i=2:6
%    measures=data.data.(fields{i});
%    split=sum(isnan(measures),2);
%    [~,I]=max(split);
%    sessionData=newmeasures(ii:5:end,:);
%    ii=ii+1;
%    data.data.(fields{i})=[data.data.(fields{i}),[sessionData(1:I-1,:);NaN(1,size(sessionData,2));sessionData(I:end,:)]];
% end

% for i=2:6
%     fields=fieldnames(data.data);  measures=data.data.(fields{i});
%     split=sum(isnan(measures),2);
%     [~,I]=max(split);
%     if i==5; measures=[measures;NaN(5,16)];end
%     Group1(:,:,i-1)=measures(1:I-1,:);
%     Group2(:,:,i-1)=measures(I+1:end,:);
% end

