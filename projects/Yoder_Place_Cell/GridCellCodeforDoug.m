% GridCellCodeforDoug
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analyses/spikeCode');
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analyses/spikeCode/read_nvt')
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/chronux_2_11/fly_track/FAnalyze/functions')
addpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')
clear;clc;close all

% PATH TO DATA
path='/Users/RyanHarvey/Downloads/2017-08-08';

% READ IN VIDEO DATA
if exist([path,filesep,'VT1.txt'],'file')==0; ReadVideoTrackingFile3(path); end
VTfile = FindFiles('VT1.txt', 'StartingDirectory',path);
tableall=textread(cell2mat(VTfile));

% LOCATE EVENTS
TSdiff=diff(tableall(:,1));
SE=find(TSdiff>std(TSdiff)*1.5);
StartofRec=[1;SE(1)+1;SE(2)+1;SE(3)+1];
EndofRec=[SE(1);SE(2);SE(3);length(tableall)];

StartofRec=tableall(StartofRec,1);
EndofRec=tableall(EndofRec,1);


for event=1:length(StartofRec)
    
    % RESTRICT DATA BY START AND END OF EVENT
    table=tableall(tableall(:,1)>StartofRec(event) & tableall(:,1)<EndofRec(event),:);
    
    
    % NEW LENGTH OF DATA            30 TIMESTAMPS = 1 SECOND
    lengthofdata=length(table);     sampleRate = 60;
    
    % DURATION OF SESSION (MIN)
    sessionduration=(lengthofdata/sampleRate)/60;
    disp(['SESSION_',num2str(event),' WAS ',num2str(sessionduration),' MIN'])
    
    % DEAL WITH NON DETECTS
    table(table==0)=NaN;
    [table,SumofNaN1,SumofNaN2,~] = NonDetects(table,lengthofdata);
    
    track_length = 61;
    
    % CHOOSE LED WITH BEST TRACKING TO VELOCITY FILTER BY
    if SumofNaN1>SumofNaN2; VelIndex=[4 5]; else VelIndex=[2 3]; end;
    
    % CALCULATE VEL TO REMOVE JUMPS IN DATA
    [vel_cmPerSec,~,~] = InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],'no',track_length,sampleRate);
    
    % find xy coords and spikes when rat is within velocity threshold ( <200)
    table = table((vel_cmPerSec<=200),:);
    
    
    filteredlengthofdata=length(table);
    disp(['VELOCITY FILTERING ', num2str(lengthofdata-filteredlengthofdata),...
        ' DATA POINTS, WHICH IS ', num2str(100-((filteredlengthofdata/lengthofdata)*100)), ' PERCENT OF YOUR DATA.']);
    
    % SMOOTHING RAW XY DATA FROM EACH LED
    padDsData=[repmat(table(1,:),30,1); table; repmat(table(end,:),30,1)]; %pad ends with last data point
    padDsData=[padDsData(:,1),runline(padDsData(:,2),10,1),runline(padDsData(:,3),10,1),runline(padDsData(:,4),10,1),runline(padDsData(:,5),10,1)];
    table=(padDsData(30+1:end-30,:));
    
    
    % GET VELOCITY
    [vel_cmPerSec,vel_abs,pixelDist] = InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],'no',track_length,sampleRate);
    
    % CONCAT TIMESTAMP/XY POSITION/ANGLE TOGETHER
    data_video_smoothfilt2=double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2)]);
    
    % FINDING ANGLE (GETS ANGLE FROM RAW XY FOR EACH LED)
    ExtractedAngle=XYangle(data_video_smoothfilt2(:,2),data_video_smoothfilt2(:,3));
    
    data_video_smoothfilt2=[data_video_smoothfilt2(2:end,:),ExtractedAngle,vel_cmPerSec]; % remove first point to fit with vel
    
    %     % EXTRACT BASIC MOVEMENT DATA
    %     BasicLoco.AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*sampleRate));
    %     BasicLoco.OverallDistTraveled=sum(vel_abs*pixelDist);
    %     BasicLoco.MeanVelocity=mean(vel_cmPerSec);
    
    tfile=FindFiles('*SS*.txt','StartingDirectory',path);
    for i=1:length(tfile)
        SpikeFile=textread(cell2mat(tfile(i)));
        SpikeFile=SpikeFile(SpikeFile(:,1)>StartofRec(event) & SpikeFile(:,1)<EndofRec(event),:);
        
        % INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
        TS = interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,1), SpikeFile, 'linear');
        X= interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,2), SpikeFile, 'linear');
        Y = interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,3), SpikeFile, 'linear');
        A = interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,4), SpikeFile, 'linear'); % *** NEED TO USE DIFFERENT INTERP FOR ANGLULAR DATA ***
        VEL = interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,5), SpikeFile, 'linear');
        
        % CONCAT AND SORT
        data_video_smoothfilt=sortrows([[TS X Y A VEL ones(size(TS,1),1)];[data_video_smoothfilt2,zeros(length(data_video_smoothfilt2),1)]],1);
        
        spks_VEL = data_video_smoothfilt(data_video_smoothfilt(:,6) == 1,:);
        
        
        [SmoothRateMap,nBinsx,nBinsy,occ,Coherence] = BinData(data_video_smoothfilt,sampleRate,spks_VEL,track_length);
        
        autocorr=SpatialAutoCorr(SmoothRateMap,length(SmoothRateMap));
        gridout=GridScore_Sinusoidal(autocorr,length(SmoothRateMap));
        
        % SAVE EVERYTHING
        Data_video_smoothfilt.(['Cell',num2str(i),'Sess',num2str(event)])=data_video_smoothfilt;
%         Spks_VEL.(['Cell',num2str(i),'Sess',num2str(event)])=spks_VEL;
        gridscore.(['Cell',num2str(i),'Sess',num2str(event)])=gridout.maxSinuGrid;
        Ratemaps.(['Cell',num2str(i),'Sess',num2str(event)])=SmoothRateMap;
        AutoCorr.(['Cell',num2str(i),'Sess',num2str(event)])=autocorr;
    end
end

% PLOT EVERYTHING
for i=1:length(tfile)
    eventi=[1,4,7,10];
    for event=1:length(StartofRec)
        data_video_smoothfilt=Data_video_smoothfilt.(['Cell',num2str(i),'Sess',num2str(event)]);
        spks_VEL = data_video_smoothfilt(data_video_smoothfilt(:,6) == 1,:);
        GS=gridscore.(['Cell',num2str(i),'Sess',num2str(event)]);
        SmoothRateMap=Ratemaps.(['Cell',num2str(i),'Sess',num2str(event)]);
        autocorr=AutoCorr.(['Cell',num2str(i),'Sess',num2str(event)]);
        
        
        fig=figure(i); fig.Color=[1 1 1]; %fig.OuterPosition=[1 6 1920 1053];
        subplot(length(StartofRec),3,eventi(event)); 
        plot(data_video_smoothfilt(:,2), data_video_smoothfilt(:,3), 'LineWidth', 1, 'color', 'k');
        hold on; axis off
        scatter(spks_VEL(:,2), spks_VEL(:,3), 10, 'filled', 'r');
        box off; axis image
        title('Spike on Path');
        
        figure (i), subplot(length(StartofRec),3,eventi(event)+1); 
        h = pcolor(SmoothRateMap);
        colormap jet; axis off; hold on; box off; set(h, 'EdgeColor', 'none'); axis image; %shading interp;
        title('Smoothed Rate Map');
        
        figure (i), subplot(length(StartofRec),3,eventi(event)+2); 
        h = pcolor(autocorr);
        colormap jet; axis off; hold on; box off; set(h, 'EdgeColor', 'none'); axis image; %shading interp;
        title(['Auto Correlation, Grid Score: ',num2str(GS)]);
    end
    p=strsplit(tfile{i},filesep);
    name=strsplit(p{end},'.'); 
    print(fig,'-dpng', '-r300',char(strcat('/Users/RyanHarvey/Downloads/2017-08-08',filesep,name{1},'GridOutput.png')))

%     print(fig,'-bestfit', '-dpdf', '-r600',char(strcat('/Users/RyanHarvey/Downloads/2017-08-08',filesep,name{1},'GridOutput.pdf')))
    close all
end





function [ SmoothRateMap,nBinsx,nBinsy,occ,Coherence] = BinData( occMatrix,sampleRate,spks_VEL,track_length)
%bindata Summary of this function goes here
%   Detailed explanation goes here

nBinsx = round(track_length/2); nBinsy = round(track_length/2);
%     mat4occ=unique(occMatrix,'rows');
MinY = min(occMatrix(:,3));
MaxY = max(occMatrix(:,3));
MinX = min(occMatrix(:,2));
MaxX = max(occMatrix(:,2));
edges{1} = linspace(MinY, MaxY, nBinsy+1);
edges{2} = linspace(MinX, MaxX, nBinsx+1);

Omatrix = hist3([occMatrix(:,3) occMatrix(:,2)],'Edges',edges);

%     Omatrix=hist3([mat4occ(:,2) mat4occ(:,3)],[nBinsy,nBinsx]);
Omatrix(end,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/sampleRate;
occ(occ<0.150)=0; % Bins with less than 150ms dropped to zero

% bin spike data
Smatrix = hist3([spks_VEL(:,3), spks_VEL(:,2)],'Edges',edges);
%     Smatrix=hist3([spks_VEL(:,2), spks_VEL(:,3)],[nBinsy nBinsx]);
Smatrix(end,:) = [];
Smatrix(:,end) = [];
% divide binned spikes by occupancy to get rate maps and store data
% (can use occ+eps instead of removing 0's)
FilledRateMatrix = Smatrix./occ;
FilledRateMatrix(isinf(FilledRateMatrix))=0;


% SMOOTH
filtWidth = [5 5]; filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout');

% COHERENCE
nonsmooth=FilledRateMatrix;
nonsmooth(isnan(nonsmooth))=0;
nonsmooth(isinf(nonsmooth))=0;

smoothed=SmoothRateMap;
smoothed(isnan(smoothed))=0;
smoothed(isinf(smoothed))=0;

Coherence=corr2(nonsmooth,smoothed);
end