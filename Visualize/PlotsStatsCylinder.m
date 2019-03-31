%% PlotsStatsCylinder
% reads in new excel data from CompileAllMatFiles and outputs stats and figures
% THIS SCRIPT SPECIFICALLY COMPILES AND RUNS STATS ON CYLINDER DATA OUTPUT FROM PAE PROJECT
%
% RYAN H 2017

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
% FIND PATH TO DATA
path='F:\Users\reharvey\Place_Cell_Data\PAE_Rat';
if ismac==1
    path='/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/PLACECELL_DATASHEETS';
end
% READ IN CYLINDER
LS17   = readtable([path,filesep,'_AllSpikeDataCylinderLS17.xlsx']);
LS19   = readtable([path,filesep,'_AllSpikeDataCylinderLS19.xlsx']);
LS21   = readtable([path,filesep,'_AllSpikeDataCylinderLS21.xlsx']);
LS23   = readtable([path,filesep,'_AllSpikeDataCylinderLS23.xlsx']);
LE2821 = readtable([path,filesep,'_AllSpikeDataCylinderLE2821.xlsx']);
LE2813 = readtable([path,filesep,'_AllSpikeDataCylinderLE2813.xlsx']);
control=[LS21;LS23;LE2821];
PAE    =[LS17;LS19;LE2813];

% EXTRACT VAR NAMES
AllVariableNames=LS17.Properties.VariableNames; AllVariableNames(1)=[];
clear LS17 LS19 LS21 LS23 LE2813 LE2821

% table to cell to get data in some sort of workable format
controlC=table2cell(control);PAEC=table2cell(PAE);

% READ IN DISPLACEMENT
LS17CD   = readtable([path,filesep,'_CylinderDisplacementLS17.xlsx']);
LS19CD   = readtable([path,filesep,'_CylinderDisplacementLS19.xlsx']);
LS21CD   = readtable([path,filesep,'_CylinderDisplacementLS21.xlsx']);
LS23CD   = readtable([path,filesep,'_CylinderDisplacementLS23.xlsx']);
LE2813CD = readtable([path,filesep,'_CylinderDisplacementLE2813.xlsx']);
LE2821CD = readtable([path,filesep,'_CylinderDisplacementLE2821.xlsx']);
controlCD=[LS21CD;LS23CD;LE2821CD];
PAECD    =[LS17CD;LS19CD;LE2813CD];
clear LS17CD LS19CD LS21CD LS23CD LE2813CD LE2821CD
controlCD=table2cell(controlCD);PAECD=table2cell(PAECD);
% % find weird name errors (future outputs shouldn't have this issue...hopefully)
% IndexC = strfind(controlC(:,1), '__');
% controlC=controlC(find(~cellfun(@isempty,IndexC)),:);
%
% IndexC = strfind(PAEC(:,1), '__');
% PAEC=PAEC(find(~cellfun(@isempty,IndexC)),:);

%% //////////////////////// ////////// Just assess first cylinder /////////////////////////////////
for s=1
    firstsessonly=1;
    if firstsessonly==1
        sessionIND=zeros(size(controlC,1),1);
        for i=1:size(controlC,1)
            temp=regexp(controlC{i,1},'\d*','Match');
            sessionIND(i,1)=str2double(temp(end));
            clear temp
        end
        controlC(:,1)=[];
        controlC=cell2mat(controlC);
        CONTROLALL.S1=controlC(sessionIND==2,:);
        
        sessionIND=zeros(size(PAEC,1),1);
        for i=1:size(PAEC,1)
            temp=regexp(PAEC{i,1},'\d*','Match');
            sessionIND(i,1)=str2double(temp(end));
            clear temp
        end
        PAEC(:,1)=[];
        PAEC=cell2mat(PAEC);
        PAEALL.S1=PAEC(sessionIND==2,:);
        
        % % % REMOVE LRATIO & ISOLATION DISTANCE - SOMETHING IS GOING ON MAKING THEM MOSTLY NANS
        CONTROLALL.S1(:,17:18)=[];
        PAEALL.S1(:,17:18)=[];
        AllVariableNames(:,17:18)=[];
        
        % Filter by Spikes
        idx=CONTROLALL.S1(:,16)>=50;
        CONTROLALL.S1=CONTROLALL.S1(idx,:);
        
        idx=PAEALL.S1(:,16)>=50;
        PAEALL.S1=PAEALL.S1(idx,:);
        
        % Filter by Peak Rate
        idx=CONTROLALL.S1(:,4)>=2;
        CONTROLALL.S1=CONTROLALL.S1(idx,:);
        
        idx=PAEALL.S1(:,4)>=2;
        PAEALL.S1=PAEALL.S1(idx,:);
        
        % Filter by Overall Firing Rate
        % idx=CONTROLALL.S1(:,5)<=10;
        % CONTROLALL.S1=CONTROLALL.S1(idx,:);
        
        % idx=PAEALL.S1(:,5)<=10;
        % PAEALL.S1=PAEALL.S1(idx,:);
        
        % Filter by field width (at least 5cm)
        % idx=CONTROLALL.S1(:,7)>=5;
        % CONTROLALL.S1=CONTROLALL.S1(idx,:);
        %
        % idx=PAEALL.S1(:,7)>=5;
        % PAEALL.S1=PAEALL.S1(idx,:);
        
        % Filter by information content
        idx=CONTROLALL.S1(:,1)>=.8;
        CONTROLALL.S1=CONTROLALL.S1(idx,:);
        
        idx=PAEALL.S1(:,1)>=.8;
        PAEALL.S1=PAEALL.S1(idx,:);
        
        
        [ AllStatsS1 ] = ScatterBox(CONTROLALL.S1,PAEALL.S1,{'Control' 'PAE'},AllVariableNames,1)
    end
end
%% //////////////////////  Only look at data that includes session 1 & 2 ////////////////////////
for s=1
    % control
    index=[];
    i=1;
    for ii=1:size(controlC,1)
        CellID1=strsplit(controlC{i,1},'__');
        CellID2=strsplit(controlC{i+1,1},'__');
        if strcmp(CellID1(1),CellID2(1))
            index=[index;1;1];
            i=i+2;
            if i>=size(controlC,1); break; end
        elseif strcmp(CellID1(1),CellID2(1))==0
            index=[index;0];
            i=i+1;
            if i>=size(controlC,1); break; end
        end
    end
    controlC=controlC(logical(index),:);
    
    % pae
    index=[];
    i=1;
    for ii=1:size(PAEC,1)
        CellID1=strsplit(PAEC{i,1},'__');
        CellID2=strsplit(PAEC{i+1,1},'__');
        if strcmp(CellID1(1),CellID2(1))
            index=[index;1;1];
            i=i+2;
            if i>=size(PAEC,1); break; end
        elseif strcmp(CellID1(1),CellID2(1))==0
            index=[index;0];
            i=i+1;
            if i>=size(PAEC,1); break; end
        end
    end
    PAEC=PAEC(logical(index),:);
    
    % GROUP DATA BY SESSION
    % control
    sessionIND=zeros(size(controlC,1),1);
    for i=1:size(controlC,1)
        temp=regexp(controlC{i,1},'\d*','Match');
        sessionIND(i,1)=str2double(temp(end));
        clear temp
    end
    Controlcluster=controlC(sessionIND==2,1);
    controlC(:,1)=[];
    controlC=cell2mat(controlC);
    CONTROLALL.S1=controlC(sessionIND==2,:);
    CONTROLALL.S2=controlC(sessionIND==3,:);
    
    % pae
    sessionIND=zeros(size(PAEC,1),1);
    for i=1:size(PAEC,1)
        temp=regexp(PAEC{i,1},'\d*','Match');
        sessionIND(i,1)=str2double(temp(end));
        clear temp
    end
    PAEcluster=PAEC(sessionIND==2,1);
    PAEC(:,1)=[];
    PAEC=cell2mat(PAEC);
    PAEALL.S1=PAEC(sessionIND==2,:);
    PAEALL.S2=PAEC(sessionIND==3,:);
    
    clear controlC controlD control PAED PAEC PAE i sessionIND
end
%% /////////////////////////////////// FILTER CELLS DOWN ///////////////////////////////////////
for s=1
    % cs1n=size(CONTROLALL.S1,1);
    % cs2n=size(CONTROLALL.S2,1);
    % ps1n=size(PAEALL.S1,1);
    % ps2n=size(PAEALL.S2,1);
    
    % LOWER AND UPPER CUT OFFS
    %Filter by Inter-Spike-Interval
    % CONTROLALL.S1=CONTROLALL.S1(CONTROLALL.S1(:,70)<2,:);
    % CONTROLALL.S2=CONTROLALL.S2(CONTROLALL.S1(:,70)<2,:);
    % PAEALL.S1=PAEALL.S1(PAEALL.S1(:,70)<2,:);
    % PAEALL.S2=PAEALL.S2(PAEALL.S1(:,70)<2,:);
    
    % Filter by Spikes
    idx=CONTROLALL.S1(:,16)>=50;
    CONTROLALL.S1=CONTROLALL.S1(idx,:);
    CONTROLALL.S2=CONTROLALL.S2(idx,:);
    Controlcluster=Controlcluster(idx,1);
    
    idx=PAEALL.S1(:,16)>=50;
    PAEALL.S1=PAEALL.S1(idx,:);
    PAEALL.S2=PAEALL.S2(idx,:);
    PAEcluster=PAEcluster(idx,1);
    
    % Filter by Peak Rate
    idx=CONTROLALL.S1(:,4)>=2;
    CONTROLALL.S1=CONTROLALL.S1(idx,:);
    CONTROLALL.S2=CONTROLALL.S2(idx,:);
    Controlcluster=Controlcluster(idx,1);
    
    idx=PAEALL.S1(:,4)>=2;
    PAEALL.S1=PAEALL.S1(idx,:);
    PAEALL.S2=PAEALL.S2(idx,:);
    PAEcluster=PAEcluster(idx,1);
    
    % Filter by Overall Firing Rate
    % idx=CONTROLALL.S1(:,5)<=10;
    % CONTROLALL.S1=CONTROLALL.S1(idx,:);
    % CONTROLALL.S2=CONTROLALL.S2(idx,:);
    % Controlcluster=Controlcluster(idx,1);
    
    % idx=PAEALL.S1(:,5)<=10;
    % PAEALL.S1=PAEALL.S1(idx,:);
    % PAEALL.S2=PAEALL.S2(idx,:);
    % PAEcluster=PAEcluster(idx,1);
    
    % % Filter by field width (at least 5cm)
    % idx=CONTROLALL.S1(:,7)>=5;
    % CONTROLALL.S1=CONTROLALL.S1(idx,:);
    % CONTROLALL.S2=CONTROLALL.S2(idx,:);
    % Controlcluster=Controlcluster(idx,1);
    %
    % idx=PAEALL.S1(:,7)>=5;
    % PAEALL.S1=PAEALL.S1(idx,:);
    % PAEALL.S2=PAEALL.S2(idx,:);
    % PAEcluster=PAEcluster(idx,1);
    
    % Filter by information content
    idx=CONTROLALL.S1(:,1)>=.8;
    CONTROLALL.S1=CONTROLALL.S1(idx,:);
    CONTROLALL.S2=CONTROLALL.S2(idx,:);
    Controlcluster=Controlcluster(idx,1);
    
    idx=PAEALL.S1(:,1)>=.8;
    PAEALL.S1=PAEALL.S1(idx,:);
    PAEALL.S2=PAEALL.S2(idx,:);
    PAEcluster=PAEcluster(idx,1);
    
    % % % REMOVE LRATIO & ISOLATION DISTANCE - SOMETHING IS GOING ON MAKING THEM MOSTLY NANS
    CONTROLALL.S1(:,17:18)=[];
    CONTROLALL.S2(:,17:18)=[];
    PAEALL.S1(:,17:18)=[];
    PAEALL.S2(:,17:18)=[];
    AllVariableNames(:,17:18)=[];
end
%% ////////////////////////////////////// Displacement ///////////////////////////////////////
% locate place cells
controlCD=cell2mat(controlCD(ismember(controlCD(:,1),Controlcluster),2:3));
PAECD=cell2mat(PAECD(ismember(PAECD(:,1),PAEcluster),2:3));
% remove low correlations
controlCD(controlCD(:,2)<.2,:)=[];
PAECD(PAECD(:,2)<.2,:)=[];
% remove NaNs
controlCD(isnan(controlCD(:,2)),:)=[];
PAECD(isnan(PAECD(:,2)),:)=[];

%%
clear ax p1 p2 Polar_Fig
Polar_Fig=figure(5);
Polar_Fig.Color=[1 1 1];
p1=polarhistogram(deg2rad(controlCD(:,1)),30,'FaceColor','k','FaceAlpha',0.8,'EdgeColor','w','Normalization','probability'); hold on
p2=polarhistogram(deg2rad(PAECD(:,1)),30,'FaceColor','red','FaceAlpha',.5,'EdgeColor','w','Normalization','probability');
ax=gca; ax.ThetaTick=[0,90,180,270]; ax.FontWeight='bold'; ax.FontSize=20; ax.RAxisLocation=45;ax.GridAlpha=.5;ax.GridColor='k';
% ax.RTick=[round(linspace(0,.16,4),2)];
% print(Polar_Fig,'-bestfit', '-dpdf', '-r300',[FigureLocation,filesep,'Polar_Fig.pdf'])
%%
% plots
figure;rose(controlCD(:,1))
figure;rose(PAECD(:,1))
% stats
circ_stats(deg2rad(controlCD(:,1)))
circ_stats(deg2rad(PAECD(:,1)))

circ_r(deg2rad(controlCD(:,1)))
circ_r(deg2rad(PAECD(:,1)))
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/FMAToolbox'))
[h,p] = ConcentrationTest([deg2rad(controlCD(:,1));deg2rad(PAECD(:,1))],[zeros(length(controlCD),1);ones(length(PAECD),1)])

circ_wwtest(deg2rad(controlCD(:,1)),deg2rad(PAECD(:,1)))
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~STATS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

% TTESTS FOR CYLINDER SESSION 1 & 2
[ AllStatsS1 ] = ScatterBox(CONTROLALL.S1,PAEALL.S1,{'Control' 'PAE'},AllVariableNames,1)

[ AllStatsS2 ] = ScatterBox(CONTROLALL.S2,PAEALL.S2,{'Control' 'PAE'},AllVariableNames,1)

%%
% [ AllStatsS1 ] = ScatterBox(CONTROLALL.S1(:,7),PAEALL.S1(:,7),{'Control' 'PAE'},AllVariableNames(7),2)
% [ AllStatsS1 ] = ScatterBox(CONTROLALL.S1(:,1),PAEALL.S1(:,1),{'Control' 'PAE'},AllVariableNames(1),2)
% [ AllStatsS1 ] = ScatterBox(CONTROLALL.S1(:,4),PAEALL.S1(:,4),{'Control' 'PAE'},AllVariableNames(4),2)
% [ AllStatsS1 ] = ScatterBox(CONTROLALL.S1(:,5),PAEALL.S1(:,5),{'Control' 'PAE'},AllVariableNames(5),2)
%
% print(figure(1), '-dpdf', '-r300',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Figures',filesep,'CyFW.pdf'])
% print(figure(2), '-dpdf', '-r300',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Figures',filesep,'CyIC.pdf'])
% print(figure(3), '-dpdf', '-r300',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Figures',filesep,'CyPR.pdf'])
% print(figure(4), '-dpdf', '-r300',['/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Figures',filesep,'CyOFR.pdf'])
%%

Polar_Fig=figure;
Polar_Fig.Color=[1 1 1];
p1=polarhistogram(deg2rad(CONTROLALL.S1(:,39)),60,'FaceColor','k','FaceAlpha',0.8,'EdgeColor','w','Normalization','probability'); hold on
p2=polarhistogram(deg2rad(PAEALL.S1(:,39)),60,'FaceColor','red','FaceAlpha',.5,'EdgeColor','w','Normalization','probability');
ax=gca; ax.ThetaTick=[0,90,180,270]; ax.FontWeight='bold'; ax.FontSize=20; ax.RAxisLocation=45;ax.GridAlpha=.5;ax.GridColor='k';
% ax.RTick=[round(linspace(0,.16,4),2)];
%%
% REPEATED MEASURES ANOVA BETWEEN STANDARD AND ROTATED SESSIONS
for i=1:length(AllVariableNames)
    control=[CONTROLALL.S1(:,i),CONTROLALL.S2(:,i)]; PAE=[PAEALL.S1(:,i),PAEALL.S2(:,i)];
    control(isnan(control(:,1)) | isnan(control(:,2)) ,:)=[]; PAE(isnan(PAE(:,1)) | isnan(PAE(:,2)) ,:)=[];
    if size(control,1)<2 || size(PAE,1)<2 ; disp([AllVariableNames(i),': Not enough Data']);continue; end
    
    groups=[repmat({'Control'},size(control,1),1);repmat({'PAE'},size(PAE,1),1)];
    
    measureVal=[control(:,1),control(:,2);PAE(:,1),PAE(:,2)];
    
    t=table(groups,measureVal(:,1),measureVal(:,2),'VariableNames',{'species','meas1','meas2'});
    Meas = table([1 2]','VariableNames',{'Measurements'});
    
    % CREATE MODELS
    % OVER ALL SESSIONS
    rm1 = fitrm(t,'meas1-meas2~species','WithinDesign',Meas);
    
    % tbl = mauchly(rm)
    % ranovatbl = ranova(rm)
    disp(strcat('--------------------------------------------------',AllVariableNames(i),'--------------------------------------------------'))
    disp('Multivariate Tests *ALL SESSIONS*')
    manovatbl=manova(rm1)
    disp('Tests of Between-Subjects Effects *ALL SESSIONS*')
    anovatbl=anova(rm1)
end

%%
Group1=cat(3,CONTROLALL.S1,CONTROLALL.S2);
Group2=cat(3,PAEALL.S1,PAEALL.S2);
for c=1:length(AllVariableNames)
    shadded_line_fig=figure; shadded_line_fig.Color=[1 1 1];shadded_line_fig.OuterPosition=[680 630 729 421];
    % CONTROL
    y=[nanmean(Group1(:,c,1)),nanmean(Group1(:,c,2))];
    SEM=[nanstd(Group1(:,c,1))/sqrt(size(Group1(:,c,1),1)),nanstd(Group1(:,c,2))/sqrt(size(Group1(:,c,2),1))];
    h1=shadedErrorBar(1:2,y,SEM,'-k',0);
    ylabel(AllVariableNames(c));
    ax1=gca; set(ax1,'XTick',[1 2],'Box','off','FontSize',20,'FontWeight','bold','LineWidth',2)
    
    hold on
    
    % TILTED
    y=[nanmean(Group2(:,c,1)),nanmean(Group2(:,c,2))];
    SEM=[nanstd(Group2(:,c,1))/sqrt(size(Group2(:,c,1),1)),nanstd(Group2(:,c,2))/sqrt(size(Group2(:,c,2),1))];
    h2=shadedErrorBar(1:2,y,SEM,'-r',1);
    set(h2.mainLine,'Color',[1 0 0])
    
    %     print(shadded_line_fig, '-dpdf', '-r300',{strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig')})
    %     print(strcat(FigureLocation,filesep,VarNames(c),'_shadded_line_fig.pdf'),'-dpdf','-r300','-bestfit')
    %         print(shadded_line_fig,'-bestfit', '-dpdf', '-r600',char(strcat('/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Figures',filesep,AllVariableNames(c),'_shadded_line.pdf')))
    cd '/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Figures'
    saveas(shadded_line_fig,char(strcat(AllVariableNames(c),'_shadded_line.jpeg')))
    close all
end

%%