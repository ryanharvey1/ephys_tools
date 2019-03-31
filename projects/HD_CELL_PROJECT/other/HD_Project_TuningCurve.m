% HD_Project_TuningCurve
% opens .r text files containing HD cell data and creates tuning curves
%
%
% Ryan E Harvey 2018
%
% cd to data and get file names
clear;clc;close all
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/CircStat2012a'))
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')

path='/Volumes/SAMSUNGUSB/Ben_HDCProject';
cd(path)
folders=dir;
folders={folders.name}';
folders(contains(folders,'._'))=[];
folders(contains(folders,'.'))=[];
folders(contains(folders,'TimeStamps'))=[];

for a=1:length(folders)
    cd(['/Volumes/SAMSUNGUSB/Ben_HDCProject/',folders{a}])
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
        
        %% CREATE TUNING CURVE
        % 6 degree bins
        da=pi/30;
        angBins=[da/2:da:2*pi-da/2];
        % Occupancy
        histAng=hist(frames(:,10),angBins);
        % Number of spikes per bin
        spkPerAng=hist(framesEXP(framesEXP(:,6)==1,10),angBins);
        % Tuning
        hdTuning=(spkPerAng./histAng)*60;
        % remove nan & inf
        hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
        
        clear da histAng spkPerAng
        
        %% COMPUTE STATS
        rlength = circ_r(angBins',hdTuning',deg2rad(6));
        [pval,z] = circ_rtest(angBins',hdTuning',deg2rad(6));
        stats(i,:)=[rlength,pval,z];
        
        clear rlength pval z
        
        %% CIRCSHIFT SO PEAK IS IN THE MIDDLE
        [M,I]=max(hdTuning);
        middlebin=round(median(1:length(hdTuning)));
        
        hdTuningshift=circshift(hdTuning,(middlebin-I)-1);
        
        %     figure;plot(hdTuning);hold on;plot([I;I],[0;M],'r')
        %     [M,I]=max(hdTuningshift);
        %     figure;plot(hdTuningshift);hold on;plot([I;I],[0;M],'r')
        %     close all
        clear M I middlebin
        
        %% Normalize & Collect tuning curve
        ALL_hdTuning(i,:)=rescale(hdTuningshift,0,1);
%                 ALL_hdTuning(i,:)=hdTuningshift;

        
    end
    %% SORT AND SAVE
    [RLength,I]=sort(stats(:,1));
    ALL_hdTuningsorted=ALL_hdTuning(I,:);
    
    currentarea=folders{a};
    index = find(isletter(currentarea), 1);
    currentarea  = currentarea(index:end);
    
    currentarea=strsplit(currentarea,'_');
    
    data.(currentarea{1}).TuningCurve=ALL_hdTuningsorted;
    data.(currentarea{1}).RLength=RLength;
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
    csvwrite(['/Volumes/SAMSUNGUSB/sortedTuningCurves',folders{a},'.csv'],data.(folders{a}).TuningCurve)
end

%% SAVE
print(tuning_Fig,'-dpng', '-r300',['/Volumes/SAMSUNGUSB',filesep,'tuning_Fig2.png'])

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

print(ALLtuning_Fig,'-dpng', '-r300',['/Volumes/SAMSUNGUSB',filesep,'ALLtuning_Fig2.png'])

%% SAVE TO CSV
csvwrite('/Volumes/SAMSUNGUSB/sortedTuningCurves2.csv',SortedTun)




