% HD_Project_TuningCurve
% opens .r text files containing HD cell data and creates tuning curves
%
%
% Ryan E Harvey 2018
%
% cd to data and get file names
clear;clc;close all

LauraDesk={'PSYCH-D385XG42'};
LabComp2={'PSYCH-DGR91CH2'};

if strcmp(getenv('computername'),LabComp2)
    addpath(genpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\CircStat2012a'))
    addpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis')
elseif strcmp(getenv('computername'),LauraDesk)
    addpath(genpath('C:\Users\Ben Clark''s Lab\Google Drive\MatlabDir\CircStat2012a'))
    addpath('C:\Users\Ben Clark''s Lab\Google Drive\MatlabDir\BClarkToolbox\Analysis')
end


if strcmp(getenv('computername'),LauraDesk)
    path='C:\Users\Ben Clark''s Lab\Google Drive\Manuscripts\In Progress\Ben_HDCProject\Data';
    cd(path)
    folders=dir;
    folders={folders.name}';
    folders(contains(folders,'._'))=[];
    folders(contains(folders,'.'))=[];
    folders(contains(folders,'TimeStamps'))=[];
elseif strcmp(getenv('computername'),LabComp2)
    path='F:\Ben_HDCProject\LM2018';
    cd(path)
    folders=dir;
    folders={folders.name}';
    folders(contains(folders,'._'))=[];
    folders(contains(folders,'.'))=[];
    folders(contains(folders,'TimeStamps'))=[];
end

for a=1:length(folders)
    cd(['C:\Users\Ben Clark''s Lab\Google Drive\Manuscripts\In Progress\Ben_HDCProject\Data\',folders{a}])
    filenames=dir('*.r');
    filenames={filenames.name}';
    filenames(contains(filenames,'._'))=[];
    
    for i=1:length(filenames)
        % OPEN TEXT FILE AND CLEAN DATA
        disp(['Running:  ',filenames{i}])
        fileID=fopen(filenames{i},'r');
        dataArray=textscan(fileID,'%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]','Delimiter',...
            '\t','TextType','string','EmptyValue',NaN,'HeaderLines',2-1,'ReturnOnError',false,'EndOfLine','\r\n');
        fclose(fileID);
        frames=[dataArray{1:end-1}];
        clear dataArray fileID
        % remove last column (all NaNs)
        frames(:,end)=[];
        
        % check firing in S2-S4
        if sum(sum(frames(:,7:9)))>0
            disp([filenames{i},' Has spikes in S2-S4'])
        end
        
        % remove non-detects in the first 4 xy columns
        frames(frames(:,2)==0 | frames(:,3)==0 | frames(:,4)==0 | frames(:,5)==0,:)=[];
        
        % expand frames to create spike binary
        spksover60hz=frames(frames(:,6)>1,:);
        
        % remove frames with > spike binary 1 for later recombining
        tempframes=frames;
        tempframes(tempframes(:,6)>1,:)=[];
        
        EXP=zeros(1,size(spksover60hz,2));
        for ii=1:size(spksover60hz,1)
            EXP=[EXP;repmat(spksover60hz(ii,:),spksover60hz(ii,6),1)];
        end
        EXP(1,:)=[];
        EXP(:,6)=ones(size(EXP,1),1);
        
        framesEXP=[tempframes;EXP];
        
        [~,I]=sort(framesEXP(:,1));
        
        framesEXP=framesEXP(I,:);
        
        
        clear tempframes EXP spksover60hz ii I
        
        %% CREATE TUNING CURVE +++++++++++++++++++++++++++++++++++++++++++
        % 6 degree bins
        da=pi/30;
        angBins=da/2:da:2*pi-da/2;
        % Occupancy
        histAng=hist(frames(:,10),angBins);
        % Number of spikes per bin
        spkPerAng=hist(framesEXP(framesEXP(:,6)==1,10),angBins);
        % Tuning
        hdTuning=(spkPerAng./histAng)*60;
        % remove nan & inf
        hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
        
        curve(i,:)=hdTuning;
        % OVERALL FIRING RATE
        OverallFR=(sum(framesEXP(:,6))/size(frames(:,6),1))*60;
        %PEAK RATE
        [peakrate,peakIdx]=max(hdTuning);
        %PREF DIR
        prefdirec=angBins(peakIdx);
%         preferred_Direction=peakIdx*6;
        %HALF WIDTH 
        [halfWidth,~]=FindHalfWidth(hdTuning);
        %DIRECTIONAL INFORMATION CONTENT
        DIC= computeDIC(histAng,hdTuning,OverallFR);
        %WITHIN QUARTER STABILITY
        [within_Coeff,within] = within_HDstability(framesEXP,frames,60,10,6);
        
        binSize=6;
        Npoints     = round(10*6/binSize); %for 6degree bins
        gw          = gausswin(Npoints,5); %alpha = 0.05
        gw          = gw/sum(gw);
        
        fig=figure; fig.Color=[1 1 1];
        for quarter=1:4
            tempTun=within.hdTuning{quarter, 1};
            %smooth curve
            l           = length(tempTun);
            hdTmp       = [tempTun tempTun tempTun];
            hdTmp       = convn(hdTmp(:),gw,'same');
            hdSmoothed  = hdTmp(l+1:2*l);
            h=plot(hdSmoothed); hold on;
            if quarter==1
                h.LineWidth=3;
                h.Color='k';
            elseif quarter ==2
                h.LineWidth=3;
                h.Color='b';
            elseif quarter ==3
                h.LineWidth=3;
                h.Color='r';
            elseif quarter==4
                h.LineWidth=3;
                h.Color=[.75 .75 .75];
            end
            box off
            set(gca,'FontWeight','Bold','FontSize',20,'LineWidth',3)
        end
        print(figure(fig), '-dpng', '-r600',['d:\Users\BClarkLab\Desktop\Laura Temp\WithinStability',filesep,['TC_within_', num2str(i),'.png']])
        close all
        
        %COHERENCE
        %Smoothing the tuning window from peyreche 
        binSize=6;
        l           = length(hdTuning);
        hdTmp       = [hdTuning hdTuning hdTuning];
        Npoints     = round(10*6/binSize); %for 6degree bins
        gw          = gausswin(Npoints,5); %alpha = 0.05
        gw          = gw/sum(gw);
        hdTmp       = convn(hdTmp(:),gw,'same');
        hdSmoothed  = hdTmp(l+1:2*l);
       
        HD_Coherence=corr2(hdSmoothed',hdTuning);
        
        %Preferred Direction Shift 
        [~,~,meanHD,N,plotIdx] = HDdrift(frames(:,10),framesEXP,frames,0);
        fig1=figure; fig1.Color=[1 1 1]; plot(meanHD,'-','Color',[.75 .75 .75],'LineWidth',3); hold on;
        y=1:length(meanHD);
        scatter(y(plotIdx==1),meanHD(plotIdx==1,:),50,'r','o','filled');
        box off
        set(gca,'FontWeight','Bold','FontSize',20,'LineWidth',3)
        print(figure(fig1), '-dpng', '-r600',['d:\Users\BClarkLab\Desktop\Laura Temp\SpikebyAngle',filesep,['TC_spikeByAngle_', num2str(i),'.png']])
        
        Measures(i,:)=[OverallFR,peak_Firing_Rate, preferred_Direction, halfWidth, DIC, within_Coeff, HD_Coherence];
        clear da histAng spkPerAng
        
        %% COMPUTE STATS
        rlength = circ_r(angBins',hdTuning',deg2rad(6));
        [pval,z] = circ_rtest(angBins',hdTuning',deg2rad(6));
        stats(i,:)=[rlength,pval,z];
        
%         Firing rate x HD polar plot for the nonsmoothed data above
        figure(i); polarplot = polar(angBins',hdTuning','b');
        set(polarplot, 'linewidth',3,'color','k'); axis off
        title(['Polor Plot, MeanVecLength: ',num2str(rlength),' Pref_Dir: ',num2str(preferred_Direction)]);
        set(0,'Showhiddenhandles','on')
        
        % ---------CODE FOR PUBLICATION FIGURE--------
        extrastuff = setdiff(get(gca,'children'),polarplot);
        delete(extrastuff)
        horizontal=line([-max(hdTuning) max(hdTuning)],[0 0]);
        vertical=line([0 0],[-max(hdTuning) max(hdTuning)]);
        set(horizontal,'linewidth',2,'color','k');
        set(vertical,'linewidth',2,'color','k');
        %---------------------------------------------
        set(figure(i),'Position',[686 325 977 619]);
        fig1 = figure(i);
%         
        clear rlength pval z
        
        %% CIRCSHIFT SO PEAK IS IN THE MIDDLE
        [M,I]=max(hdTuning);
        middlebin=round(median(1:length(hdTuning)));
        
        hdTuningshift=circshift(hdTuning,(middlebin-I)-1);
        
            figure;plot(hdTuning);hold on;plot([I;I],[0;M],'r')
            [M,I]=max(hdTuningshift);
            figure;plot(hdTuningshift);hold on;plot([I;I],[0;M],'r')
            close all
        clear M I middlebin
        
        %% Normalize & Collect tuning curve
        ALL_hdTuning(i,:)=rescale(hdTuningshift,0,1);
        ALL_hdTuning(i,:)=hdTuningshift;

        
    end
    %% SORT AND SAVE
    [RLength,I]=sort(stats(:,1));
    ALL_hdTuningsorted=ALL_hdTuning(I,:);
    ALL_measuressorted=Measures(I,:);
    
    currentarea=folders{a};
    index = find(isletter(currentarea), 1);
    currentarea  = currentarea(index:end);
    
    currentarea=strsplit(currentarea,'_');
    data.(currentarea{1}).varnames={'OverallFR','peak_Firing_Rate', 'preferred_Direction', 'halfWidth', 'DIC', 'within_Coeff','HD_Coherence'};
    data.(currentarea{1}).TuningCurve=curve;
    data.(currentarea{1}).RLength=stats;
    data.(currentarea{1}).Measures=Measures;
end

%% GRAPH Tuning curves sorted by RLength arranged by brain area
folders=fieldnames(data);

tuning_Fig=figure; tuning_Fig.Color=[1 1 1];

for a=1:length(folders)
    subplot(2,2,a)
    imagesc(data.(folders{a}).TuningCurve);
    axis xy; colormap jet; box off;
    ylabel('Cells')
    set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
    title(folders{a})
    csvwrite(['F:\Ben_HDCProject\LM2018\taubecontrols',folders{a},'.csv'],data.(folders{a}).TuningCurve)
end

%% SAVE
print(tuning_Fig,'-dpng', '-r300',['F:\Ben_HDCProject\LM2018\',filesep,'tuning_Fig2.png'])

%% GRAPH All tuning curves in one matrix sorted by RLength

ALL_R=[];
ALL_TUN=[zeros(1,size(data.(folders{1}).TuningCurve,2))];
for a=1:length(folders)
    ALL_R=[ALL_R;data.(folders{a}).RLength];
    ALL_TUN=[ALL_TUN;data.(folders{a}).TuningCurve];
end
ALL_TUN(1,:)=[];

[SortedR,I]=sort(ALL_R);
SortedTun=ALL_TUN(I,:);

ALLtuning_Fig=figure; ALLtuning_Fig.Color=[1 1 1];

imagesc(SortedTun);
axis xy; colormap jet; box off;
ylabel('Cells')
set(gca,'XTick',[1,30.7,60],'XTickLabel',[-180 0 180],'FontSize',20,'FontWeight','bold')
title('ATN,MEC,PaS,PoS')

print(ALLtuning_Fig,'-dpng', '-r300',['F:\Ben_HDCProject\LM2018',filesep,'ALLtuning_Fig2.png'])

%% SAVE TO CSV
% csvwrite('F:\Ben_HDCProject\LM2018\sortedTuningCurves2.csv',SortedTun)




