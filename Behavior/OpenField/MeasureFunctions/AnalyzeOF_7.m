
%%
%%%%%%%%%%% Plot Paths with home base outline and cue location
load('params_V8');
param_idx=params.PathName;
params.cueCenter{1}=[];

c=[params.transcue{19}(4,:);params.transcue{19}(2,:);params.transcue{19}(1,:);params.transcue{19}(3,:)];
cue = polyshape(c);
[cueX,cueY]=centroid(cue);
     

for j=1:length(params.transcoords)
    
    if  isempty(params.rateMap{j})
        continue
    end
    
    figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k','LineWidth',2);
    
    hold on; plot(params.transcoords{j}(:,1),params.transcoords{j}(:,2),'--k','LineWidth',2)
    
    for hb=1:size(params.HBcoords{j},1)
        hold on; plot(params.HBcoords{j}{hb,1},params.HBcoords{j}{hb,2},'r','LineWidth',2)
    end
    
    
    if contains(param_idx{j},'day2_lgOF')
        hold on; plot(cueX,cueY,'*r', 'MarkerSize',25,'LineWidth',3) 
    else
        hold on; o=plot(cue);
        o.FaceColor=[0 0 0];
        o.FaceAlpha=1;
    end
    
    axis image
    axis off
    title(erase(param_idx{j},["D:\Maria\OF_Excel\",".xlsx","_"]))
    
    print(gcf,'-dpng', '-r300',['D:\Maria\Figures\',erase(param_idx{j},["D:\Maria\OF_Excel\",".xlsx"]),'.png'])
    close all
    
end



%%
% %%%%%%%%%%%%% Ratemaps

for i =1:length(params.transcoords)

    if isempty(params.rateMap{i})
        continue
    end

map=params.rateMap{i};

figure;
PerfectCircRateMap(map,1);
 title(erase(param_idx{i},["D:\Maria\OF_Excel\",".xlsx","_"]))

print(gcf,'-dpng', '-r300',['D:\Maria\Figures\',erase(param_idx{i},["D:\Maria\OF_Excel\",".xlsx"]),'heatmaps.png'])
    close all

end


%%
%%Plot cue and circle
figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k','LineWidth',2);
 c=[params.transcue{19}(4,:);params.transcue{19}(2,:);params.transcue{19}(1,:);params.transcue{19}(3,:)];

  cue = polyshape(c);
  hold on; cue=plot(cue)
  cue.FaceColor=[0 0 0];
  cue.FaceAlpha=1;

hold on; plot(params.transcue{19}(:,1),params.transcue{19}(:,2))
axis image
axis off

%%Plot Circle Without Cue (red Star)

figure; plot(cos(0:2*pi/1000:2*pi)*101,sin(0:2*pi/1000:2*pi)*101,'-k','LineWidth',2);
  hold on; plot(cueX,cueY,'*r', 'MarkerSize',25,'LineWidth',3) 
axis image
axis off

%% Path Length 

pathL_tg1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
pathL_wt1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

pathL_tg2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
pathL_wt2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});


% RUN MIXED ANOVA
RL_anova([pathL_tg1 pathL_tg2],[pathL_wt1 pathL_wt2],'Path Length')


%% BEESWARM PLOTS FOR PL, NUM STOPS, NUM RUNS, & PATH CIRCUITY 
%%Idx the two groups for the two days
NumRuns_tg1=vertcat(params.NumRuns{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
NumRuns_wt1=vertcat(params.NumRuns{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

NumRuns_tg2=vertcat(params.NumRuns{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
NumRuns_wt2=vertcat(params.NumRuns{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});

RL_anova([NumRuns_wt1 NumRuns_wt2],[NumRuns_tg1 NumRuns_tg2],'Num Runs')

medianNumRuns_tg1=median(NumRuns_tg1);
medianNumRuns_tg2=median(NumRuns_tg2);
medianNumRuns_wt1=median(NumRuns_wt1);
medianNumRuns_wt2=median(NumRuns_wt2);

%%Varargins for the function and to customize the graph
data=[NumRuns_tg1; NumRuns_wt1;NumRuns_tg2;NumRuns_wt2];
dayIdx = ['Day1';'Day2']; %%labels x axis
xValues_idx=[1;2];%the values of the x axis
distBeesIdx=[ones(1,(24)) (ones(1,(24))+1)]';
yLabel_idx= 'Number of Runs';
catIdx=[ones(1,12) (ones(1,12)+1) ones(1,12) (ones(1,12)+1)]';

%%Creates a graph for Tg v Wt of each group 
f=figure; f.Color=[1 1 1];
plotSpread(data,'categoryIdx',catIdx,'xNames',dayIdx,'xValues',xValues_idx,...
    'categoryColors',{'r',[.25,.25,.25]},'DistributionIdx',distBeesIdx,'yLabel',yLabel_idx,...
    'distributionMarkers','o');
 hold on;
 plot([medianNumRuns_tg1,medianNumRuns_tg2], 'LineWidth',2,'Color','r')
 hold on;
 plot([medianNumRuns_wt1,medianNumRuns_wt2], 'LineWidth',2,'Color',[.25 .25 .25])
 
 %%%Editing Figure
legend('Tg','Wt','Tg Median','Wt Median')
title('Number of Runs','fontweight','bold')
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',24)

%%Idx the two groups for the two days
numStops_tg1=vertcat(params.NumStops{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
numStops_wt1=vertcat(params.NumStops{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

numStops_tg2=vertcat(params.NumStops{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
numStops_wt2=vertcat(params.NumStops{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});

RL_anova([numStops_tg1 numStops_tg2],[numStops_wt1 numStops_wt2],'Number of Stops')

mediannumStops_tg1=median(numStops_tg1);
mediancnumStops_tg2=median(numStops_tg2);
mediannumStops_wt1=median(numStops_wt1);
mediannumStops_wt2=median(numStops_wt2);

%%Varargins for the function and to customize the graph
data=[numStops_tg1; numStops_wt1;numStops_tg2;numStops_wt2];
dayIdx = ['Day1';'Day2']; %%labels x axis
xValues_idx=[1;2];%the values of the x axis
distBeesIdx=[ones(1,(24)) (ones(1,(24))+1)]';
yLabel_idx= 'Number of Stops';
catIdx=[ones(1,12) (ones(1,12)+1) ones(1,12) (ones(1,12)+1)]';

%%Creates a graph for Tg v Wt of each group 
f=figure; f.Color=[1 1 1];
plotSpread(data,'categoryIdx',catIdx,'xNames',dayIdx,'xValues',xValues_idx,...
    'categoryColors',{'r',[.25,.25,.25]},'DistributionIdx',distBeesIdx,'yLabel',yLabel_idx,...
    'distributionMarkers','o');
 hold on;
 plot([mediannumStops_tg1,mediannumStops_tg2], 'LineWidth',2,'Color','r')
 hold on;
 plot([mediannumStops_wt1,mediannumStops_wt2], 'LineWidth',2,'Color',[.25 .25 .25])
 
 %%%Editing Figure
legend('Tg','Wt','Tg Median','Wt Median')
title('Circuity','fontweight','bold')
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',24)

%%Idx the two groups for the two days
circuity_tg1=horzcat(params.meanCirc{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')})';
circuity_wt1=horzcat(params.meanCirc{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')})';

circuity_tg2=horzcat(params.meanCirc{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')})';
circuity_wt2=horzcat(params.meanCirc{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')})';

RL_anova([circuity_tg1 circuity_tg2],[circuity_wt1 circuity_wt2],'Segment Circuity')

mediancircuity_tg1=median(circuity_tg1);
mediancircuity_tg2=median(circuity_tg2);
mediancircuity_wt1=median(circuity_wt1);
mediancircuity_wt2=median(circuity_wt2);

%%Varargins for the function and to customize the graph
data=[circuity_tg1; circuity_wt1;circuity_tg2;circuity_wt2];
dayIdx = ['Day1';'Day2']; %%labels x axis
xValues_idx=[1;2];%the values of the x axis
distBeesIdx=[ones(1,(24)) (ones(1,(24))+1)]';
yLabel_idx= ['Circuity'];
catIdx=[ones(1,12) (ones(1,12)+1) ones(1,12) (ones(1,12)+1)]';


%%Creates a graph for Tg v Wt of each group 
f=figure; 
axesf=gca
h=plotSpread(data,'categoryIdx',catIdx,'xNames',dayIdx,'xValues',xValues_idx,...
    'categoryColors',{'r',[.25,.25,.25]},'DistributionIdx',distBeesIdx,'yLabel',yLabel_idx,...
    'distributionMarkers','o');
 hold on;
 plot([mediancircuity_tg1,mediancircuity_tg2], 'LineWidth',2,'Color','r')
 hold on;
 plot([mediancircuity_wt1,mediancircuity_wt2], 'LineWidth',2,'Color',[.25 .25 .25])
 
 %%%Editing Figure
legend('Tg','Wt','Tg Median','Wt Median')
title('Circuity','fontweight','bold')
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',24)

f.Color=[1 1 1];

%%

%Create example of dwell times with cue added

[zones] = createZones([0,0], 202,'fig',1,'numquad',16);
hold on; 
plot(cueX,cueY,'*r', 'MarkerSize',25,'LineWidth',3) 
hold on; 
plot(cueX*-1,cueY*-1,'*b', 'MarkerSize',25,'LineWidth',3) 

%%DWELL TIME
dwellQuad_tg1=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
dwellQuad_wt1=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

dwellQuad_tg2=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
dwellQuad_wt2=vertcat(params.dwellQuad{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});

%Histograms for quadrant dwell time
bincenter=[0:360/16:360]-22.5;
bincenter=bincenter(2:end);

%Histograms for quadrant dwell time
f=figure; f.Color=[1 1 1];
subplot(1,2,1)
bar(sum(dwellQuad_wt1,1)/sum(dwellQuad_wt1(:)),'FaceColor','k','FaceAlpha',.5);
% title('Day 1 "Cue" Quadrant Dwell Time for TG (red) and WT (gray)')
ylim([0 .6])
xlabel('Quadrant')
ylabel('Normalized Dwell Time')
 hold on 
 bar(sum(dwellQuad_tg1,1)/sum(dwellQuad_tg1(:)),'FaceColor','r','FaceAlpha',.5);
 set(gca,'FontSize',24,'FontWeight','bold','FontName','Helvetica','LineWidth',2,'XTick',1:2:16)%'XTick',1:16,'XTickLabel',round(bincenter),
box off 

subplot(1,2,2)
bar(sum(dwellQuad_wt2,1)/sum(dwellQuad_wt2(:)),'FaceColor','k','FaceAlpha',.5);
% title('Day 2 "No Cue" Quadrant Dwell Time for TG (red) and WT (gray)')
ylim([0 .6])
xlabel('Quadrant')
ylabel('Normalized Dwell Time')
 hold on 
 bar(sum(dwellQuad_tg2,1)/sum(dwellQuad_tg2(:)),'FaceColor','r','FaceAlpha',.5);
  set(gca,'FontSize',24,'FontWeight','bold','FontName','Helvetica','LineWidth',2,'XTick',1:2:16); %'XTick',1:16,'XTickLabel',round(bincenter),
box off 

%%

%%%Bar Chart for Scenarios
f=figure; f.Color=[1 1 1];
bar([1,2,3,4,5],[5,2,5,0,0],'k','FaceAlpha',.5)
hold on;
bar([1,2,3,4,5],[1,4,3,1,3],'r','FaceAlpha',.5)

ylabel('Number of Rats')
legend('Wt','Tg')
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',24,'YTick',0:1:5)

box off

%% Segment PATH Length CDF plot


segPL_tg1=horzcat(params.segPL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1'),:});
segPL_wt1=horzcat(params.segPL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1'),:});

segPL_tg2=horzcat(params.segPL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2'),:});
segPL_wt2=horzcat(params.segPL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2'),:});

f=figure; f.Color=[1 1 1]; 
CoolHistogram(cell2mat(segPL_tg1)',cell2mat(segPL_wt1)',25,'Segment Path Length - Day 1')
ylim([0 .7])

[ AllStats ] = CDFplots(cell2mat(segPL_wt1)',cell2mat(segPL_tg1)',{'WT','TG'},{'Segment Length (cm)'},2)
% 
% f=figure; f.Color=[1 1 1]; 
% CoolHistogram(cell2mat(segPL_tg2)',cell2mat(segPL_wt2)',25,'Segment Path Length - Day 2')
% ylim([0 .7])
hold on; 
[ AllStats ] = CDFplots(cell2mat(segPL_wt2)',cell2mat(segPL_tg2)',{'WT','TG'},{'Segment Length (cm)'},2)

%% Compute number of HBs and # of stops in hb

for j=1:size(param_idx,1)
    
    if  isempty(params.rateMap{j})
        continue
    end
    
   params.numHB{j}=size(params.HBcoords{j},1);

end
%% PROXIMITY CDF 

prox_tg1=horzcat(params.proximity{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1'),:});
prox_wt1=horzcat(params.proximity{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1'),:});

prox_tg2=horzcat(params.proximity{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2'),:});
prox_wt2=horzcat(params.proximity{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2'),:});
% 
% f=figure; f.Color=[1 1 1]; 
% CoolHistogram(cell2mat(prox_tg1)',cell2mat(prox_wt1)',25,'Segment Path Length - Day 1')
% ylim([0 .7])

[ AllStats ] = CDFplots(prox_wt1',prox_tg1',{'WT','TG'},{'Stop Proximity to HB (cm) - Day 1'},2)
% 
% f=figure; f.Color=[1 1 1]; 
% CoolHistogram(cell2mat(segPL_tg2)',cell2mat(segPL_wt2)',25,'Segment Path Length - Day 2')
% ylim([0 .7])

[ AllStats ] = CDFplots(prox_wt2',prox_tg2',{'WT','TG'},{'Stop Proximity to HB (cm) - Day 2'},2)

%% 

numHB_tg1=vertcat(params.numHB{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
numHB_wt1=vertcat(params.numHB{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

numHB_tg2=vertcat(params.numHB{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
numHB_wt2=vertcat(params.numHB{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});


for i=1:5
    numHB_wt1_cat(i)=sum(numHB_wt1==i);
end

for i=1:5
    numHB_tg1_cat(i)=sum(numHB_tg1==i);
end

for i=1:5
    numHB_wt2_cat(i)=sum(numHB_wt2==i);
end

for i=1:5
    numHB_tg2_cat(i)=sum(numHB_tg2==i);
end

f=figure; f.Color=[1 1 1]; 
subaxis(2,1,1, 'Spacing', 0.2, 'Padding', 0);
barh(1:5,numHB_wt1_cat,'FaceColor','k','FaceAlpha',.5)
hold on
barh(1:5,numHB_tg1_cat,'FaceColor','r','FaceAlpha',.5)
box off
ylim([0 5.5])
xlabel('# Animals')
ylabel('# Home Bases')
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)
title('Day 1')

subaxis(2,1,2, 'Spacing', 0.2, 'Padding', 0);
barh(1:5,numHB_wt2_cat,'FaceColor','k','FaceAlpha',.5)
hold on
barh(1:5,numHB_tg2_cat,'FaceColor','r','FaceAlpha',.5)
box off
ylim([0 5.5])
xlabel('# Animals')
ylabel('# Home Bases')
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)
title('Day 2')




%%
stopsNearHB_tg1=vertcat(params.stopsNearHB{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
stopsNearHB_wt1=vertcat(params.stopsNearHB{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});

stopsNearHB_tg2=vertcat(params.stopsNearHB{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
stopsNearHB_wt2=vertcat(params.stopsNearHB{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});

f=figure; f.Color=[1 1 1]; 
scatter(stopsNearHB_wt1,stopsNearHB_wt2,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
hold on
scatter(stopsNearHB_tg1,stopsNearHB_tg2,'r','Filled');


RL_anova([stopsNearHB_wt1,stopsNearHB_wt2],[stopsNearHB_tg1,stopsNearHB_tg2],'Number of Stops')

mstopsNearHB_tg1=mean(stopsNearHB_tg1);
mstopsNearHB_tg2=mean(stopsNearHB_tg2);
mstopsNearHB_wt1=mean(stopsNearHB_wt1);
mstopsNearHB_wt2=mean(stopsNearHB_wt2);

%%Varargins for the function and to customize the graph
data=[stopsNearHB_tg1; stopsNearHB_wt1;stopsNearHB_tg2;stopsNearHB_wt2];
dayIdx = ['Day1';'Day2']; %%labels x axis
xValues_idx=[1;2];%the values of the x axis
distBeesIdx=[ones(1,(24)) (ones(1,(24))+1)]';
yLabel_idx= 'Number of Stops';
catIdx=[ones(1,12) (ones(1,12)+1) ones(1,12) (ones(1,12)+1)]';

%%Creates a graph for Tg v Wt of each group 
f=figure; f.Color=[1 1 1];
plotSpread(data,'categoryIdx',catIdx,'xNames',dayIdx,'xValues',xValues_idx,...
    'categoryColors',{'r',[.25,.25,.25]},'DistributionIdx',distBeesIdx,'yLabel',yLabel_idx,...
    'distributionMarkers','o');
 hold on;
 plot([mstopsNearHB_tg1,mstopsNearHB_tg2], 'LineWidth',3,'Color','r')
 hold on;
 plot([mstopsNearHB_wt1,mstopsNearHB_wt2], 'LineWidth',3,'Color',[.25 .25 .25])
 
 %%%Editing Figure
legend('Tg','Wt','Tg Median','Wt Median')
title('Number of Stops in Home Base','fontweight','bold')
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',24)

%% 


Scenario=[1 2 2 2 2 3 3 3 4 5 5 5 1 1 1 1 1 2 2 3 3 3 3 3];
Group=[zeros(1,12) ones(1,12)];


[TgWtChi,chi2,p,labels] = crosstab(Scenario,Group) %between Genotype


% f=figure; f.Color=[1 1 1]; 
% CoolHistogram(numHB_wt1,numHB_tg1,5,'# of Home Bases - Day 1')
% set(gca,'XTick',1:3)
% 
% f=figure; f.Color=[1 1 1]; 
% CoolHistogram(numHB_wt2,numHB_tg2,5,'# of Home Bases - Day 2')
% 











%%%%%%%%%%%%%%CODE GRAVEYARD%%%%%%%%%%%%%%%%%%
% % % % % %%"BEE SWARM" PLOT of Path Length
% % % % %
% % % % % param_idx=params.PathName;
% % % % % addpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\Visualize\plotSpread')
% % % % %
% % % % % %%Creates variables of each group and each day
% % % % %     pathL_tg1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day1')});
% % % % %     pathL_wt1=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day1')});
% % % % %
% % % % %     pathL_tg2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Tg') & contains(param_idx,'day2')});
% % % % %     pathL_wt2=vertcat(params.pathL{contains(param_idx,'lgOF') & contains(param_idx,'Wt') & contains(param_idx,'day2')});
% % % % %
% % % % % %%Varargins for the function and to customize the graph
% % % % %         %%Important Note: shownMW creates a line for either mean median or both, '1' creates for
% % % % %         %%both(used in line 116)
% % % % % data=[pathL_tg1; pathL_wt1;pathL_tg2;pathL_wt2];
% % % % % group=[ones(1,12) (ones(1,12)+1)]';
% % % % % catIdx=[group; group];
% % % % % dayIdx = ['Day1';'Day2']; %%labels x axis
% % % % % xValues_idx=[1;2];%the values of the x axis
% % % % % distBeesIdx=[zeros(24,1); ones(24,1)];
% % % % % yLabel_idx= ['Path Length Score'];
% % % % %
% % % % %
% % % % % %%Plotting
% % % % % h=figure;
% % % % % h=plotSpread(data,'categoryIdx',catIdx,'xNames',dayIdx,'xValues',xValues_idx,...
% % % % %     'categoryColors',{'g',[.5,.5,.5]},'DistributionIdx',distBeesIdx,'yLabel',yLabel_idx,'showMM',0,...
% % % % %     'distributionMarkers','.');
% % % % % set(h,'MarkerSize',25)
% % % % % title('Path Length')
% % % % % legend('Tg','Wt')
% % % % %
% % % % % box off
% % % % %


