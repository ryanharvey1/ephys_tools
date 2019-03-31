%% PlotsStats
% reads in new excel data from CompileAllMatFiles and outputs stats and figures
clc, clear, close all

addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/MatlabStatsUofG'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\MatlabStatsUofG'));
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/RC_notBoxPlot'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\RC_notBoxPlot'));

addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\CircStat2012a'));
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/CircStat2012a'));
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/images'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\images'));
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis'));


% ###################
placecells=true;
interneurons=false;
% ###################

path='F:\Users\reharvey\Place_Cell_Data\PAE_Rat';
if ismac==1
    path='/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/PLACECELL_DATASHEETS';
end
tableexample = readtable([path,filesep,'_AllSpikeDataRH11.xlsx']);
tableexample(:,1)=[];
RH11 = xlsread([path,filesep,'_AllSpikeDataRH11.xlsx']);
RH16 = xlsread([path,filesep,'_AllSpikeDataRH16.xlsx']);
RH13 = xlsread([path,filesep,'_AllSpikeDataRH13.xlsx']);
RH14 = xlsread([path,filesep,'_AllSpikeDataRH14.xlsx']);
LS17 = xlsread([path,filesep,'_AllSpikeDataLS17.xlsx']);
LS19 = xlsread([path,filesep,'_AllSpikeDataLS19.xlsx']);
LS21 = xlsread([path,filesep,'_AllSpikeDataLS21.xlsx']);
LS23 = xlsread([path,filesep,'_AllSpikeDataLS23.xlsx']);
LE2821 = xlsread([path,filesep,'_AllSpikeDataLE2821.xlsx']);
LE2813 = xlsread([path,filesep,'_AllSpikeDataLE2813.xlsx']);


control=[RH13;RH14;LS21;LS23;LE2821]; %PhControl=control;
PAE=[RH11;RH16;LS17;LS19;LE2813]; %PhPAE=PAE;
numofcells=length(control)+length(PAE);
clear RH11 RH16 RH13 RH14 LS17 LS19 LS21 LS23 LE2821 LE2813

% % EXTRACT BASIC MOVEMENT
% BasicMovControl=control(:,end-2:end); % * change if new vars added
% BasicMovPAE=PAE(:,end-2:end);
% % remove NaNs (you will only have one set of values per session )
% BasicMovControl(isnan(BasicMovControl(:,1)),:)=[];
% BasicMovPAE(isnan(BasicMovPAE(:,1)),:)=[];
% % delete movement vars from main data sets (for the above reason)
% control(:,end-2:end)=[]; % * change if new vars added
% PAE(:,end-2:end)=[];


%% SEP CELL TYPE
for i=1
if placecells==true
    
%     %Filter by PH precession
%     control=control(control(:,68)<0,:);
%     PAE=PAE(PAE(:,68)<0,:);
    
%     %Filter by Inter-Spike-Interval
%     control=control(control(:,76)<2,:);
%     PAE=PAE(PAE(:,76)<2,:);
    
    % Filter by Spikes
    control=control(control(:,16)>=50,:);
    PAE=PAE(PAE(:,16)>=50,:);
    
    % Filter by Peak Rate
    control=control(control(:,4)>=2,:);
    PAE=PAE(PAE(:,4)>=2,:);
    
    %    Filter by mean firing rate
    %     control=control(control(:,5)>=.25,:);
    %     PAE=PAE(PAE(:,5)>=.25,:);
    
%     %     % Filter by Overall Firing Rate
%     control=control(control(:,5)<=10,:);
%     PAE=PAE(PAE(:,5)<=10,:);
    
    %     % Filter by field distance to wall (get rid of fields at ends of track)
    %     control=control(control(:,6)>4,:);
    %     PAE=PAE(PAE(:,6)>4,:);
    
%     % Filter by field width (at least 2 bins)
%     control=control(control(:,7)>=2,:);
%     PAE=PAE(PAE(:,7)>=2,:);
    
    % Filter by coherence
%     control=control(control(:,2)>=0.7,:);
%     PAE=PAE(PAE(:,2)>=0.7,:);
    
    % Filter by information content
    control=control(control(:,1)>=0.8,:);
    PAE=PAE(PAE(:,1)>=0.8,:);
end
if interneurons==true
    % Spikes
    control=control(control(:,16)>=50,:);
    PAE=PAE(PAE(:,16)>=50,:);
    
    % Info Content
    control=control(control(:,1)<=0.4808,:);
    PAE=PAE(PAE(:,1)<=0.4808,:);
    
    % Filter by Peak Rate
    control=control(control(:,4)>=1,:);
    PAE=PAE(PAE(:,4)>=1,:);
    
    % Average firing rate
    control=control(control(:,5)>=7,:);
    PAE=PAE(PAE(:,5)>=7,:);
end
% Convert directionality index to positive values
control(:,22)=abs(control(:,22));
PAE(:,22)=abs(PAE(:,22));

% % CUMULATIVE FREQ OF directionality index
% subplot()
% [f,x] = ecdf(control(:,22));
% p4=plot(x,f);
% set(p4,'LineWidth',4,'Color','k')
% hold on
% [f,x] = ecdf(PAE(:,22));
% p4_2=plot(x,f);
% set(p4_2,'LineWidth',4,'Color','r')
% box off
% legend([p4 p4_2],'Control','PAE','FontSize',12,'Location','best')
% xlabel('Directionality Index');ylabel('Cumulative Frequency')
% 
% % CUMULATIVE FREQ OF displacement
% [f,x] = ecdf(control(:,77));
% p4=plot(x,f);
% set(p4,'LineWidth',4,'Color','k')
% hold on
% [f,x] = ecdf(PAE(:,77));
% p4_2=plot(x,f);
% set(p4_2,'LineWidth',4,'Color','r')
% box off
% legend([p4 p4_2],'Control','PAE','FontSize',12,'Location','best')
% xlabel('Displacement');ylabel('Cumulative Frequency')
% 
% % CUMULATIVE FREQ OF directional correlation
% [f,x] = ecdf(control(:,77));
% p4=plot(x,f);
% set(p4,'LineWidth',4,'Color','k')
% hold on
% [f,x] = ecdf(PAE(:,77));
% p4_2=plot(x,f);
% set(p4_2,'LineWidth',4,'Color','r')
% box off
% legend([p4 p4_2],'Control','PAE','FontSize',12,'Location','best')
% xlabel('Displacement');ylabel('Cumulative Frequency')


disp(['Analyzing: ',num2str(length(control)+length(PAE)),' out of ',num2str(numofcells),' cells.']);

AllVariableNames=tableexample.Properties.VariableNames;
[ AllStats ] = ScatterBox(control,PAE,{'Control' 'PAE'},AllVariableNames,1)

end
%% Phase Locking
% Filter by Spikes
% PhControl=PhControl(PhControl(:,16)>=50,:); PhPAE=PhPAE(PhPAE(:,16)>=50,:);
% % Filter by Peak Rate
% PhControl=PhControl(PhControl(:,4)>=.5,:); PhPAE=PhPAE(PhPAE(:,4)>=.5,:);
% Calculate stats and create plots
% [PhaseStats] = PhaseLock(PhControl(:,(41:46)),PhPAE(:,(41:46)),PhControl(:,(29:34)),PhPAE(:,(29:34)));
[PhaseStats] = PhaseLock(control(:,(41:46)),PAE(:,(41:46)),control(:,(29:34)),PAE(:,(29:34)));
%%
% WORKING*******
% Frequency Histograms showing the distribution of phase precession
% Conw=control(:,68); Conw(isnan(Conw))=[]; Conw=zscore(Conw);
% PAEw=PAE(:,68); PAEw(isnan(PAEw))=[]; PAEw=zscore(PAEw);
% figure(5); subplot(1,2,1); h1=histogram(control(:,68));
% xlim([-1 1]); title('Control'); ylabel('Number of Fields')
% subplot(1,2,2); h2=histogram(PAE(:,68));
% xlim([-1 1]);title('PAE');ylabel('Number of Fields');
% set(h1,'FaceColor',[0.247 0.247 0.247])
% set(h2,'FaceColor',[0.247 0.247 0.247])
% 
% clear data;data.Control=control(:,69); data.PAE=PAE(:,69);
% [freq]=ECDF_plot(data,AllVariableNames(69));
close all

figure;
subplot(1,2,1)
histogram(control(:,23));hold on; histogram(control(control(:,29)<.05,23));ylim([0 50])
subplot(1,2,2)
histogram(PAE(:,23));hold on; histogram(PAE(PAE(:,29)<.05,23));ylim([0 50])

[ AllStats ] = ScatterBox(control(control(:,29)<.05,23),PAE(PAE(:,29)<.05,23),{'Control' 'PAE'},AllVariableNames(23),2)

clear data;data.Control=control(control(:,29)<.05,23); data.PAE=PAE(PAE(:,29)<.05,23);
[freq]=ECDF_plot(data,AllVariableNames(23));

Polar_Fig=figure;
Polar_Fig.Color=[1 1 1];
p1=polarhistogram(deg2rad(control(control(:,29)<.05,41)),30,'FaceColor','k','FaceAlpha',0.8,'EdgeColor','w','Normalization','probability'); hold on
p2=polarhistogram(deg2rad(PAE(PAE(:,29)<.05,41)),30,'FaceColor','red','FaceAlpha',.5,'EdgeColor','w','Normalization','probability');
ax=gca; ax.ThetaTick=[0,90,180,270]; ax.FontWeight='bold'; ax.FontSize=20; ax.RAxisLocation=45;ax.GridAlpha=.5;ax.GridColor='k';
% ax.RTick=[round(linspace(0,.16,4),2)];
%%
M=median([control(:,77);PAE(:,77)]);

controlabove=control(control(:,77)>M,:);
controlbelow=control(control(:,77)<M,:);
PAEabove=PAE(PAE(:,77)>M,:);
PAEbelow=PAE(PAE(:,77)<M,:);
[ AllStatscontrol ] = ScatterBox(controlabove,controlbelow,{'Control' 'PAE'},AllVariableNames,1)
[ AllStatscontrol ] = ScatterBox(PAEabove,PAEbelow,{'Control' 'PAE'},AllVariableNames,1)


freqR=find(contains(AllVariableNames,'DepthofModulation'));
% freqP=find(contains(AllVariableNames,'Rayleigh_Theta'));

% controlabove=controlabove(controlabove(:,freqP)<.05,freqR);
% controlbelow=controlbelow(controlbelow(:,freqP)<.05,freqR);
% PAEabove=PAEabove(PAEabove(:,freqP)<.05,freqR);
% PAEbelow=PAEbelow(PAEbelow(:,freqP)<.05,freqR);

[ AllStatscontrol ] = ScatterBox(controlabove(:,freqR),PAEabove(:,freqR),{'Control' 'PAE'},AllVariableNames(freqR),2)
[ AllStatsPAE ] = ScatterBox(controlbelow(:,freqR),PAEbelow(:,freqR),{'Control' 'PAE'},AllVariableNames(freqR),2)
print(figure(3),'-bestfit', '-dpdf', '-r300',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Figures',filesep,'DepofModBelow.pdf'])


%% Directionality index
clear data;data.Control=control(:,22); data.PAE=PAE(:,22);
[freq]=ECDF_plot(data,AllVariableNames(22));
[ AllStatscontrol ] = ScatterBox(data.Control,data.PAE,{'Control' 'PAE'},AllVariableNames(22),2)

%%
% % REMOVE PVALUES for Graphs and stats
% AllVariableNames(:,29:34)=[]; % remove pvalues
% control(:,29:34)=[];
% PAE(:,29:34)=[];

% [ AllStats ] = ScatterBox(control,PAE,{'Control' 'PAE'},AllVariableNames,1)

% for i=[7,1,4,5]
% % [ AllStats ] = ScatterBox(control(:,i),PAE(:,i),{'Control' 'PAE'},AllVariableNames(i),2)
% end

% print(figure(10),'-bestfit', '-dpdf', '-r300',['/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Figures',filesep,'OverallFR.pdf'])
% belowc=control(control(:,77)<16,:);
% abovec=control(control(:,77)>16,:);
% 
% belowp=PAE(PAE(:,77)<16,:);
% abovep=PAE(PAE(:,77)>16,:);
% 
% [ AllStats ] = ScatterBox(abovep,belowp,{'above' 'below'},AllVariableNames,1)
% [ AllStats ] = ScatterBox(abovec,belowc,{'above' 'below'},AllVariableNames,1)
% 
% [ AllStats ] = ScatterBox(belowc,belowp,{'belowc' 'belowp'},AllVariableNames,1)
% [ AllStats ] = ScatterBox(abovec,abovep,{'abovec' 'abovep'},AllVariableNames,1)
% 


% CORRELATIONS
% correlated to running speed
% CorrelationsControl=[];
% for i=1:size(control,2)
%     Correlation=corr2(control(:,i),control(:,61));
%     CorrelationsControl=[CorrelationsControl;AllVariableNames(i),Correlation];
% end
% CorrelationsPAE=[];
% for i=1:size(PAE,2)
%     Correlation=corr2(PAE(:,i),PAE(:,61));
%     CorrelationsPAE=[CorrelationsPAE;AllVariableNames(i),Correlation];
% end


%% WHERE OLD CODE GOES TO DIE \/\/\/
%     % correlation coefficient
%     data = nan(max(size(Group1DegWork,1),size(Group2DegWork,1)),2);
%     data((1:size(Group1DegWork,1)),1)=Group1DegWork; data((1:size(Group2DegWork,1)),2)=Group2DegWork;
%
%     [rho pval]=circ_corrcc(deg2rad(data(:,1)),deg2rad(data(:,2)));
%

% close all
% % Simple Stats
% variables={'Information Content','Sparsity','Peak Rate','Number of Spikes','DirectionalityIndex','Dis From Track End',...
%     'RLengthDelta','RLengthTheta','RLengthAlpha','RLengthBeta','RLengthGamma','RLengthHighGamma',...
%     'RayleighZDelta','RayleighZTheta','RayleighZAlpha','RayleighZBeta','RayleighZGamma','RayleighZHighGamma',...
%     'MeanPhaseDelta','MeanPhaseTheta','MeanPhaseAlpha','MeanPhaseBeta','MeanPhaseGamma','MeanPhaseHighGamma'};

% % Position of above variables in matrix
% variable_number=[1,3,4,16,22,6,...
%     23,24,25,26,27,28,...
%     35,36,37,38,39,40,...
%     41,42,43,44,45,46];
% %how big is my screen?
% scrsz = get(groot,'ScreenSize');
% fig1=figure('Position',[1.5, 1.5, scrsz(3), scrsz(4)]);
%
% for i=1:length(variables)
%     fig1; subplot(4,6,i)
%     bar([nanmean(control(:,variable_number(i))),nanmean(PAE(:,variable_number(i)))],'FaceColor',[.7 .7 .7]);
%     hold on
%     % plot(1,control(:,variable_number(i)),'.','MarkerFaceColor','k'); plot(2,PAE(:,variable_number(i)),'.','MarkerFaceColor','k');
%     if i==7 || i==8 || i==9 || i==10 || i==11 ||i==12; ylim([0 .15]);% max(max([control(:,22:27);PAE(:,22:27)]))])
%     elseif i==13 || i==14 || i==15 || i==16 || i==17 ||i==18; ylim([0  4]);%max(max([control(:,34:39);PAE(:,34:39)]))]);
%     elseif i==19 || i==20 || i==21 || i==22 || i==23 ||i==24 ; ylim([0  200]);%max(max([control(:,40:45);PAE(:,40:45)]))]),
%     end
%     hold on
%     errorbar(1:2,[nanmean(control(:,variable_number(i))),nanmean(PAE(:,variable_number(i)))],...
%         [((std(control(:,variable_number(i)))/(sqrt(length(control(:,variable_number(i))))))) ,...
%         ((std(PAE(:,variable_number(i)))/(sqrt(length(PAE(:,variable_number(i)))))))],'.','Color','k');
%     title(variables(i),'FontSize',14);
%     xlabel('Control     PAE','FontSize',14);
%     ylabel(variables(i),'FontSize',14);
%     set(gca,'box','off')
%     hold off
%
%     [~,~,ci,stats] = ttest2(control(:,variable_number(i)),PAE(:,variable_number(i)));
%     [p,~,statz] = ranksum(control(:,variable_number(i)),PAE(:,variable_number(i)));
%     staaaats(i,:)=([strcat(cellstr(variables(i)),': t(',num2str(stats.df),')=',num2str(statz.zval),', p=',(num2str(p))),ci']);
% end

% figRaster=figure(3);
% for i=40:45
% %     figRaster; subplot(6,1,i)
%     raster(control(:,i))
% end


%%

% ii=1;
% iii=7;
% for i=40:45
% [PValControl(ii),ZValueControl(ii)] = circ_rtest(deg2rad(control(:,i)));
% [PValPAE(ii),ZValuePAE(ii)] = circ_rtest(deg2rad(PAE(:,i)));
%
% MeanResultantPhaseControl(ii)=circ_r(deg2rad(control(:,i))); % R LENGTH
% MeanResultantPhasePAE(ii)=circ_r(deg2rad(PAE(:,i))); % R LENGTH
%
% ii=ii+1;
%
% figure(2); subplot(2,6,i)
% rose(control(:,i))
% subplot(2,6,iii)
% rose(PAE(:,i))
% iii=iii+1;
% end
% controlstats=[PValControl;ZValueControl;MeanResultantPhaseControl]
% PAEstats=[PValPAE;ZValuePAE;MeanResultantPhasePAE]
%
% FigRose=figure(2);
% ii=1;
% iii=7;
% variables={'DELTA','THETA','ALPHA','BETA','GAMMA','HI-GAMMA'};
% for i=40:45
%
%     FigRose; subplot(2,6,ii)
%     rose(control(:,i))
%     title([variables(ii),' R: ',num2cell(circ_r(deg2rad(control(:,i))))],'FontSize',9);
%
%     FigRose; subplot(2,6,iii)
%     rose(PAE(:,i))
%     title([variables(ii),' R: ',num2cell(circ_r(deg2rad(PAE(:,i))))],'FontSize',9);
%
%     ii=ii+1;
%     iii=iii+1;
% end
%
% FigHist=figure(3);
% ii=1;
% variables={'DELTA','THETA','ALPHA','BETA','GAMMA','HI-GAMMA'};
% for i=40:45
%     FigHist; subplot(2,3,ii)
%     h1=histogram(control(:,i),25);
%     hold on
%     h1.Normalization = 'countdensity';
%     hold on
%     h2=histogram(PAE(:,i),25);
%     hold on
%     h2.Normalization = 'countdensity';
%     xlim([50 300])
%     ylim([0 16])
%     title(variables(ii),'FontSize',12);
%     xlabel('Degrees')
%     ylabel('Normalized Spikes Per Degree','FontSize',12)
%     ii=ii+1;
% end
%
%
% % saveas(fig1,[path filesep '_PlaceCellPlots.tiff']);
% disp('DONE')
%
%
% %%
% % SPIKES ON PHASE!!!
% close all
% variables={'DELTA','THETA','ALPHA','BETA','GAMMA','HI-GAMMA'};
% for iii=1:2
%     if iii==1;group=control; elseif iii==2;group=PAE; end;
%     ii=1;
%     for i=41:46
%         if iii==1; figure(i); elseif iii==2; figure(i+6); end;
%         % plot wave
%         steps=(2*pi)/(length(group)-1);
%         t = 0:steps:(2*pi);
%         y = sin(t);
%         hh=plot(t,y,'k'); hold on
%         hh.LineWidth=5;
%         set(gca, 'XTick', []);
%         set(gca, 'YTick', []);
%         box off
%         title(variables(ii),'FontSize',12); ii=ii+1;
%         % scatter spikes
%         spikes=sin(deg2rad(group(:,i))');
%         h=scatter(deg2rad(group(:,i))',spikes,'Marker','x');
%         h.MarkerEdgeColor=[1 0 0];
%         h.MarkerFaceColor=[1 0 0];
%         h.LineWidth = 8;
%         axis([0 2*pi -1.2 1.2])
%     end
% end

%


% % Phase Locking Z Scores
% figure(2)
% bar([nanmean(control(:,33)),nanmean(PAE(:,33)); nanmean(control(:,34)),nanmean(PAE(:,34));...
%     nanmean(control(:,35)),nanmean(PAE(:,35)); nanmean(control(:,36)),nanmean(PAE(:,36));...
%     nanmean(control(:,37)),nanmean(PAE(:,37)); nanmean(control(:,38)),nanmean(PAE(:,38))]);
% hold on; ylabel('Phase Locking (Z Scores)','FontSize',14);
%
% % Phase Locking R Length
% figure(3)
% bar([nanmean(control(:,21)),nanmean(PAE(:,21)); nanmean(control(:,22)),nanmean(PAE(:,22));...
%     nanmean(control(:,23)),nanmean(PAE(:,23)); nanmean(control(:,24)),nanmean(PAE(:,24));...
%     nanmean(control(:,25)),nanmean(PAE(:,25)); nanmean(control(:,26)),nanmean(PAE(:,26))]);
% hold on; ylabel('Phase Locking (R Length)','FontSize',14);
%
% % Phase Locking Degree
% figure(4)
% bar([nanmean(control(:,39)),nanmean(PAE(:,39)); nanmean(control(:,40)),nanmean(PAE(:,40));...
%     nanmean(control(:,41)),nanmean(PAE(:,41)); nanmean(control(:,42)),nanmean(PAE(:,42));...
%     nanmean(control(:,43)),nanmean(PAE(:,43)); nanmean(control(:,44)),nanmean(PAE(:,44))]);
% hold on; ylabel('Phase Locking (Degree)','FontSize',14);



% saveas(fig1,[path '_BasicFilteredStats.tiff']);
% %%
% clear i
% %######################################################################################################
%
% % FILTER OUT <50 SPIKES
% control=control(control(:,15)>=50,:);
% PAE=PAE(PAE(:,15)>=50,:);
% fig1;
% scrsz = get(groot,'ScreenSize');
% fig2=figure('Position',[1, 281, scrsz(3), scrsz(4)/4]);
% for i=1:length(variables)
% fig2; subplot(1,6,i)
% bar([nanmean(control(:,variable_number(i))),nanmean(PAE(:,variable_number(i)))],'FaceColor',[.7 .7 .7]);
% if i==1, ylim ([0 1.5]), elseif i==2, ylim ([0 .90]), elseif i==3, ylim ([0 2.5]), elseif i==4, ylim ([0 1700]),...
% elseif i==5, ylim ([0 .5]), elseif i==6, ylim ([0 .25]), end
% hold on
% errorbar(1:2,[nanmean(control(:,variable_number(i))),nanmean(PAE(:,variable_number(i)))],...
%     [((std(control(:,variable_number(i)))/(sqrt(length(control(:,variable_number(i))))))) ,...
%     ((std(PAE(:,variable_number(i)))/(sqrt(length(PAE(:,variable_number(i)))))))],'.','Color','k');
% title(variables(i),'FontSize',14);
% xlabel('Control              PAE','FontSize',14);
% ylabel(variables(i),'FontSize',14);
% set(gca,'box','off')
% end
% clear i
% %######################################################################################################
%
% % FILTER OUT <.5 Information Content
% control=control(control(:,1)>=.50,:);
% PAE=PAE(PAE(:,1)>=.50,:);
% fig1;
% fig2;
% scrsz = get(groot,'ScreenSize');
% fig3=figure('Position',[1, 2, scrsz(3), scrsz(4)/4]);
% for i=1:length(variables)
% fig3; subplot(1,6,i)
% bar([nanmean(control(:,variable_number(i))),nanmean(PAE(:,variable_number(i)))],'FaceColor',[.7 .7 .7]);
% if i==1, ylim ([0 1.5]), elseif i==2, ylim ([0 .90]), elseif i==3, ylim ([0 2.5]), elseif i==4, ylim ([0 1700]),...
% elseif i==5, ylim ([0 .5]), elseif i==6, ylim ([0 .25]), end
% hold on
% errorbar(1:2,[nanmean(control(:,variable_number(i))),nanmean(PAE(:,variable_number(i)))],...
%     [((std(control(:,variable_number(i)))/(sqrt(length(control(:,variable_number(i))))))) ,...
%     ((std(PAE(:,variable_number(i)))/(sqrt(length(PAE(:,variable_number(i)))))))],'.','Color','k');
% title(variables(i),'FontSize',14);
% xlabel('Control              PAE','FontSize',14);
% ylabel(variables(i),'FontSize',14);
% set(gca,'box','off')
% end
% clear i
%
% % FILTER OUT ><.5 MEAN VECTOR LENGTH
% VectBelow5 = filtspikes(filtspikes(:,18)<=.5,:);
% VectAbove5 = filtspikes(filtspikes(:,18)>=.5,:);
%
% VectAbove5ICAbove5=VectAbove5(VectAbove5(:,1)>=.5,:);
% meanMVL=mean(VectAbove5ICAbove5(:,18))
% meanIC=mean(VectAbove5ICAbove5(:,1))
% meanDirIC=mean(VectAbove5ICAbove5(:,20))
%
%  %%
% % FILTER OUT ><5HZ
% HiPassfilt = filtspikes(filtspikes(:,4)<=5,:);
% LowPassfilt = filtspikes(filtspikes(:,4)>=5,:); % interneurons
%
% % can you get the mean and median for information content, sparsity, & firing rate
% fprintf('Mean for information content: %d. median for information content: %d below 5hz\n',nanmean(HiPassfilt(:,1)),nanmedian(HiPassfilt(:,1)));
% fprintf('Mean for information content: %d. median for information content: %d above 5hz\n',nanmean(LowPassfilt(:,1)),nanmedian(LowPassfilt(:,1)));
% fprintf('Mean for sparsity: %d. median for sparsity: %d below 5hz\n',nanmean(HiPassfilt(:,2)),nanmedian(HiPassfilt(:,2)));
% fprintf('Mean for sparsity: %d. median for sparsity: %d above 5hz\n',nanmean(LowPassfilt(:,2)),nanmedian(LowPassfilt(:,2)));
% fprintf('Mean for firing rate: %d. Median for firing rate: %d. below 5hz\n',nanmean(HiPassfilt(:,4)),nanmedian(HiPassfilt(:,4)));
% fprintf('Mean for firing rate: %d. Median for firing rate: %d. above 5hz\n',nanmean(LowPassfilt(:,4)),nanmedian(LowPassfilt(:,4)));
%
% %  % of below 5Hz cells that have a greater than 0.7 information content?
% percent=(length(HiPassfilt(HiPassfilt(:,1)>=0.7,:))/length(HiPassfilt(:,1)))*100;
% fprintf('Percent of below 5Hz cell that have greater than 0.7 info content %d\n',percent);
%
% fprintf('Mean for mean vector length: %d. Median for mean vector length: %d.\n',nanmean(VectAbove5(:,18)),nanmedian(VectAbove5(:,18)));
% percent=length(VectAbove5)/length(filtspikes)*100;
% fprintf('Percent of cells above mean vector length of .5 is: %d\n',percent);
%
% figure(3)
% hist(filtspikes(:,18),100);
% xlabel('mean vector length');
%
%
% %%
% % PLOTS
% % firing rate
% figure(i),subplot(3,2,1)
%     hist(HiPassfilt(:,4),100);
%     xlabel('Firing Rate Below 5hz');
%     xlim([0 5])
%     ylim([0 350]);
% figure(i),subplot(3,2,2)
%     hist(LowPassfilt(:,4),100);
%     xlabel('Firing Rate Above 5hz');
%     xlim([5 100])
%     ylim([0 350]);
%
% % information content
% figure(i),subplot(3,2,3)
%     hist(HiPassfilt(:,1),100);
%     xlabel('information content Below 5hz');
%     xlim([0 3]);
%     ylim([0 550]);
%
% figure(i),subplot(3,2,4)
%     hist(LowPassfilt(:,1),100);
%     xlabel('information content Above 5hz ');
%     xlim([0 3]);
%     ylim([0 550]);
% % sparsity
% figure(i),subplot(3,2,5)
%     hist(HiPassfilt(:,2),100);
%     xlabel('Sparsity Below 5hz');
%     xlim([0 1]);
%     ylim([0 100]);
% figure(i),subplot(3,2,6)
%     hist(LowPassfilt(:,2),100);
%     xlabel('Sparsity Above 5hz');
%     xlim([0 1]);
%     ylim([0 100]);
%
% % cd('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/All_Mat_output_10_14_16')
%
%
