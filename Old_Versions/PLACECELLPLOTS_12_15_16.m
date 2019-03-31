%% PLACECELLPLOTS_12_15_16 
% reads in new excel data from CompileAllMatFiles and outputs stats and figures
clc, clear, close all  
tableexample = readtable('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH11.xlsx');
tableexample(:,1)=[];

placecells=true;
interneurons=false;
% ###################
path='F:\Users\reharvey\Place_Cell_Data\PAE_Rat';
    RH11 = xlsread('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH11.xlsx');
    RH16 = xlsread('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH16.xlsx');
    RH13 = xlsread('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH13.xlsx');
    RH14 = xlsread('F:\Users\reharvey\Place_Cell_Data\PAE_Rat\_AllSpikeDataRH14.xlsx');
    
    control=[RH13;RH14];
    PAE=[RH11;RH16];
    %
    if placecells==true
        % Filter by Spikes
        control=control(control(:,16)>=50,:);
        PAE=PAE(PAE(:,16)>=50,:);
        
        % Filter by Peak Rate
        control=control(control(:,4)>=1,:);
        PAE=PAE(PAE(:,4)>=1,:);
        
        % Filter by information content
        control=control(control(:,1)>=0.454702890198576,:);
        PAE=PAE(PAE(:,1)>=0.454702890198576,:);
    end
    if interneurons==true
%         Filter for interneuron
        control=control(control(:,1)<=.20,:);
        PAE=PAE(PAE(:,1)<=.20,:);
        
        control=control(control(:,16)>=50,:);
        PAE=PAE(PAE(:,16)>=50,:);
        
        control=control(control(:,5)>=7,:);
        PAE=PAE(PAE(:,5)>=7,:);
    end
    
    control(:,22)=abs(control(:,22));
    PAE(:,22)=abs(PAE(:,22));
    
    clear RH11 RH16 RH13 RH14
    disp(['Amount of cells currently analyzed:  ',num2str(length(control)+length(PAE))]);
%######################################################################################################

close all
% Simple Stats
variables={'Information Content','Sparsity','Peak Rate','Number of Spikes','DirectionalityIndex','Dis From Track End',...
    'RLengthDelta','RLengthTheta','RLengthAlpha','RLengthBeta','RLengthGamma','RLengthHighGamma',...
    'RayleighZDelta','RayleighZTheta','RayleighZAlpha','RayleighZBeta','RayleighZGamma','RayleighZHighGamma',...
    'MeanPhaseDelta','MeanPhaseTheta','MeanPhaseAlpha','MeanPhaseBeta','MeanPhaseGamma','MeanPhaseHighGamma'};

% Position of above variables in matrix
variable_number=[1,3,4,16,22,6,... 
                23,24,25,26,27,28,...
                35,36,37,38,39,40,...
                41,42,43,44,45,46];
%how big is my screen?           
scrsz = get(groot,'ScreenSize'); 
fig1=figure('Position',[1.5, 1.5, scrsz(3), scrsz(4)]);

for i=1:length(variables)
fig1; subplot(4,6,i)
bar([nanmean(control(:,variable_number(i))),nanmean(PAE(:,variable_number(i)))],'FaceColor',[.7 .7 .7]);
hold on
% plot(1,control(:,variable_number(i)),'.','MarkerFaceColor','k'); plot(2,PAE(:,variable_number(i)),'.','MarkerFaceColor','k');
if i==7 || i==8 || i==9 || i==10 || i==11 ||i==12; ylim([0 .15]);% max(max([control(:,22:27);PAE(:,22:27)]))])
elseif i==13 || i==14 || i==15 || i==16 || i==17 ||i==18; ylim([0  4]);%max(max([control(:,34:39);PAE(:,34:39)]))]);
elseif i==19 || i==20 || i==21 || i==22 || i==23 ||i==24 ; ylim([0  200]);%max(max([control(:,40:45);PAE(:,40:45)]))]),
end
hold on
errorbar(1:2,[nanmean(control(:,variable_number(i))),nanmean(PAE(:,variable_number(i)))],...
    [((std(control(:,variable_number(i)))/(sqrt(length(control(:,variable_number(i))))))) ,...
    ((std(PAE(:,variable_number(i)))/(sqrt(length(PAE(:,variable_number(i)))))))],'.','Color','k');
title(variables(i),'FontSize',14);
xlabel('Control     PAE','FontSize',14);
ylabel(variables(i),'FontSize',14);
set(gca,'box','off')
hold off

[~,~,ci,stats] = ttest2(control(:,variable_number(i)),PAE(:,variable_number(i)));
[p,~,statz] = ranksum(control(:,variable_number(i)),PAE(:,variable_number(i)));
staaaats(i,:)=([strcat(cellstr(variables(i)),': t(',num2str(stats.df),')=',num2str(statz.zval),', p=',(num2str(p))),ci']);
end

% figRaster=figure(3);
% for i=40:45
% %     figRaster; subplot(6,1,i)
%     raster(control(:,i))
% end


%%
addpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\CircStat2012a')

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
FigRose=figure(2);
ii=1;
iii=7;
variables={'DELTA','THETA','ALPHA','BETA','GAMMA','HI-GAMMA'};
for i=40:45
    
FigRose; subplot(2,6,ii)
rose(control(:,i))
title([variables(ii),' R: ',num2cell(circ_r(deg2rad(control(:,i))))],'FontSize',9);

FigRose; subplot(2,6,iii)
rose(PAE(:,i))
title([variables(ii),' R: ',num2cell(circ_r(deg2rad(PAE(:,i))))],'FontSize',9);

ii=ii+1;
iii=iii+1;
end

FigHist=figure(3);
ii=1;
variables={'DELTA','THETA','ALPHA','BETA','GAMMA','HI-GAMMA'};
for i=40:45
FigHist; subplot(2,3,ii)
h1=histogram(control(:,i),25);
hold on
h1.Normalization = 'countdensity';
hold on
h2=histogram(PAE(:,i),25);
hold on
h2.Normalization = 'countdensity';
xlim([50 300])
ylim([0 16])
title(variables(ii),'FontSize',12);
xlabel('Degrees')
ylabel('Normalized Spikes Per Degree','FontSize',12)
ii=ii+1;
end


% saveas(fig1,[path filesep '_PlaceCellPlots.tiff']);
disp('DONE')


%%
% SPIKES ON PHASE!!!    
close all
variables={'DELTA','THETA','ALPHA','BETA','GAMMA','HI-GAMMA'};
for iii=1:2
    if iii==1;group=control; elseif iii==2;group=PAE; end;
    ii=1;
    for i=41:46
        if iii==1; figure(i); elseif iii==2; figure(i+6); end;
        % plot wave
        steps=(2*pi)/(length(group)-1);
        t = 0:steps:(2*pi);
        y = sin(t);
        hh=plot(t,y,'k'); hold on
        hh.LineWidth=5;
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        box off
        title(variables(ii),'FontSize',12); ii=ii+1;
        % scatter spikes
        spikes=sin(deg2rad(group(:,i))');
        h=scatter(deg2rad(group(:,i))',spikes,'Marker','x');
        h.MarkerEdgeColor=[1 0 0];
        h.MarkerFaceColor=[1 0 0];
        h.LineWidth = 8;
        axis([0 2*pi -1.2 1.2])
    end
end

%%


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
%%
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
