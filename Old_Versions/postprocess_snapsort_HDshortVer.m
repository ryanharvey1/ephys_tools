
%%%%%%%% DATA EXTRACTION & POSTPROCESSING from Neuralynx SNAP sort file %%%%%%%%%%%%%
%
% THIS SCRIPT CREATES FIGURES AND STATS FOR PLACE & HD CELL DATA
% CAN HANDLE CIRCULAR ARENA & LINEAR TRACK DATA
%
% INPUTS: PATH TO RAW NON-PROCESSED DATA (NON-P FOLDER) / TRACK LENGTH / MAZE-TYPE / FIGURES
%   -path: path that leads to subject raw data folder (aka F:\P30_Recordings\LB01\2016-08-02_14-16-08_SS1_1)
%   -track_Length: longest width of maze (cm)
%   -linear_track: 'yes' or 'no'
%   -figures: 0 or 1
%
% OUTPUTS: RATE MAP / SPIKE ON PATH / TUNING CURVE / POLAR PLOT / R VS. L SPIKES / STATS / LFP ANALYSIS
% Ryan Harvey: 12/21/17
% EDITS By Laura Berkowitz: 02/11/2018

%    __    __    __    __    __    __    __    __    __    __    __
% __/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__
%   \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/


function postprocess_snapsort_HDshortVer(path,track_length,linear_track)
close all;clc
% ADD TOOLBOXES TO PATH
com=pwd;com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'BClarkToolbox',filesep,'Analyses',filesep,'spikeCode',filesep,'MEX'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analyses',filesep,'spikeCode'],...
    [basedir,filesep,'chronux_2_11'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'LFP'],...
    [basedir,filesep,'CircStat2012a'],...
    [basedir,filesep,'FMAToolbox'],...
    genpath([basedir,filesep,'neuralynxmex']),...
    [basedir,filesep,'BClarkToolbox',filesep,'Analyses',filesep,'spikeCode',filesep,'read_nvt']);

% sets path to snap folder & looks within for mat files
cd([path,filesep,'SNAPSorterResults']);
file=struct2table(dir( '**/*.mat'));
t=table2cell(file(:,1));
file=t(~contains(t,'_info.mat'),1);
info=t(contains(t,'_info.mat'),1);

%'LRatio','Tightness','ShortISI','Incompleteness','IsolationDistance','NumberOfSpikes','StationInTime','TempMatch','BDistanceClust','BDistanceSpike'
ns=1;grade=[];clusterquality=[];
for m=1:length(file)
    load(file{m});load(info{m})
    grade=[grade;final_grades'];
    clusterquality=[clusterquality;grades(:,1:10)];
    numspk=unique(output(:,2));
    for n=1:length(numspk)
        ID{ns,1}=orig_filename;
        avgwave{ns,1}=means{1,n};
        S{ns,1}=output(output(:,2)==n,1); ns=ns+1;
    end
end
disp(path)
if ~exist('S','var')
    disp('No Clusters... Skipping')
    return
end
disp(['    ',num2str(length(S)),' Cells'])

% CALCULATE WAVEFORM PROPERTIES
prop=waveformprop(avgwave); %RYANS MEASURES

% READ VIDEO FILE
% Creates video tracking text file if one doesn't exist
if exist([path,filesep,'VT1.txt'],'file')==0; ReadVideoTrackingFile3(path); end

% locates newly created VT1.txt file
fileID=fopen([path,filesep,'VT1.txt'],'r');
dataArray=textscan(fileID,'%f%f%f%f%f%[^\n\r]','Delimiter','\t','TextType','string','EmptyValue',NaN,'ReturnOnError',false);
fclose(fileID);
video=[dataArray{1:end-1}];
clear dataArray fileID

% SPLIT UP SESSIONS BY EVENT FILES
[ timestamps,StartofRec,EndofRec ]=EventSplit(path);if isempty(timestamps);clear timestamps;end

% SET UP VARIABLE NAMES & DATA STRUCTURE

varnames={'bins_Sampled','mean_vector_length','half_width','Directional_IC','peak_Firing_Rate','Overall_Firing_Rate','NSpikes','preferred_Direction',...
    'LRatio','Tightness','ShortISI','Incompleteness','IsolationDistance',...
    'NumberOfSpikes','StationInTime','TempMatch','BDistanceClust','BDistanceSpike',...
    'WaveformPeak','WaveformValley','spikewidth','halfwidth','peak2valleyratio',...
    'peak2valleyslope','SpikeAmplitude','rLength_95th','within_Stability','Coherence',...
    'FieldWidth','infieldFR','outfieldFR','InformationContent','borderScore','Field2Wall',...
    'PeakRate','sparsity','NumbActiveBins','thetaindex','thetapeak','MeanThetaFreq','MeanOverallPow','ThetaRatio',...
    'phaselockR','phaselockPval','HD_Coherence','PDDMax','PDDMin'};

ratID=strsplit(path,filesep);
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).varnames=varnames;
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).Spikes=cellfun(@(x) x*1000000,S,'un',0); %cell array of timestamps
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).avgwave=avgwave; %from snap

% ADD ENTIRE VIDEO TO DATA
[table,~,~,~]=NonDetects(video);
% fix tracker jumps
[table(:,2),table(:,3)]=trackerjumps(table(:,2),table(:,3));
[table(:,4),table(:,5)]=trackerjumps(table(:,4),table(:,5));
padDsData=[repmat(table(1,:),30,1); table; repmat(table(end,:),30,1)];
padDsData=[padDsData(:,1),smooth(padDsData(:,2),10),smooth(padDsData(:,3),10),smooth(padDsData(:,4),10),smooth(padDsData(:,5),10)];
table=(padDsData(30+1:end-30,:));
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).frames...
    =double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2)]); %put complete processed video data into matrix

% ADD EVENT TO DATA
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).events=[StartofRec;EndofRec]; %StartofRec= [timeStart;timeEnd timeStart;timeEnd] aka for 2 sessions

% START OF EVENT LOOP 1:nSESSIONS
for event=1:length(StartofRec)
    
    % RESTRICT DATA BY START AND END OF EVENT
    if exist('timestamps','var')==1
        try
            table=video(video(:,1)>StartofRec(event) & video(:,1)<EndofRec(event),:);
        catch
            disp(['Event File Failure in Session',num2str(event)]);
            continue
        end
    else
        table=video;
    end
    
    % NEW LENGTH OF DATA
    lengthofdata=length(table);sampleRate=30;
    
    % DURATION OF SESSION (MIN)
    sessionduration=(lengthofdata/sampleRate)/60;
    sessionduration_sec=(lengthofdata/sampleRate); %For Shuffling 
    disp(['SESSION_',num2str(event),' WAS ',num2str(sessionduration),' MIN'])
    if sessionduration<2;disp('SESSION TOO SHORT...SKIPPING');continue;end
    
    % DEAL WITH NON DETECTS
    [table,SumofNaN1,SumofNaN2,~]=NonDetects(table);
    
    % CHOOSE LED WITH BEST TRACKING TO VELOCITY FILTER BY
    if SumofNaN1>SumofNaN2; VelIndex=[4 5]; else; VelIndex=[2 3]; end
    
    % fix tracker jumps
    [table(:,2),table(:,3)]=trackerjumps(table(:,2),table(:,3));
    [table(:,4),table(:,5)]=trackerjumps(table(:,4),table(:,5));
     
    % SMOOTHING RAW XY DATA FROM EACH LED
    padDsData=[repmat(table(1,:),30,1); table; repmat(table(end,:),30,1)]; %pad ends with last data point
    padDsData=[padDsData(:,1),smooth(padDsData(:,2),10),smooth(padDsData(:,3),10),smooth(padDsData(:,4),10),smooth(padDsData(:,5),10)];
    table=(padDsData(30+1:end-30,:)); %remove pad
    
    % GET VELOCITY
    [vel_cmPerSec,vel_abs,pixelDist]=InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,sampleRate);
    
    % CONCAT TIMESTAMP/XY POSITION/ANGLE TOGETHER
    data_video=double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2)]);
    
    % GET HEAD ANGLE BASED OFF LEDs
    ExtractedAngle=XYangleLED(table(:,2),table(:,3),table(:,4),table(:,5));
    data_video=[data_video,ExtractedAngle,[vel_cmPerSec; vel_cmPerSec(end)]]; % remove first point to fit with vel
    
    % EXTRACT BASIC MOVEMENT DATA
    BasicLoco.AverageAnglePerSec=rad2deg(circ_mean((abs(diff(ExtractedAngle)))*sampleRate));
    BasicLoco.OverallDistTraveled=sum(vel_abs*pixelDist);
    BasicLoco.MeanVelocity=mean(vel_cmPerSec);
    data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).BasicLoco=BasicLoco;
    
    % START OF MAIN SPIKE TO PATH AND DATA ANALYSIS LOOP
    for i=1:length(S)
        cell=strsplit(ID{i},filesep);
        disp([cell{end},'   Cell: ',num2str(i)])
        
        SpikeFile=S{i}*1000000;
        
        % RESTRICT DATA BY START AND END OF EVENT
        if exist('timestamps','var')==1
            SpikeFile=SpikeFile(SpikeFile(:,1)>StartofRec(event) & SpikeFile(:,1)<EndofRec(event),:);
        end
        
        % FIND INDEX FOR VELOCITY FILTER (<2cm/second for 0 second)
        
        in=data_video(:,5)<0; %changed from <2s 11-FEb-2018
        dsig=diff([0 (abs(in')>=eps) 0]);
        startIndex=find(dsig>0);endIndex=find(dsig<0)-1;
        stringIndex=(endIndex-startIndex+1>=60);
        startIndex=startIndex(stringIndex);
        endIndex=endIndex(stringIndex);
        indices=zeros(1,max(endIndex)+1);
        indices(startIndex)=1;indices(endIndex+1)=indices(endIndex+1)-1;
        indices=find(cumsum(indices));
        in=zeros(length(in),1);in(indices',1)=1;
        data_video_nospk=[data_video,in];
        
        % INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
        X=interp1(data_video_nospk(:,1),data_video_nospk(:,2),SpikeFile,'linear');
        Y=interp1(data_video_nospk(:,1),data_video_nospk(:,3),SpikeFile,'linear');
        A=interp1(data_video_nospk(:,1),data_video_nospk(:,4),SpikeFile,'linear');
        VEL=interp1(data_video_nospk(:,1),data_video_nospk(:,5),SpikeFile,'linear');
        VELidx=interp1(data_video_nospk(:,1),data_video_nospk(:,6),SpikeFile,'linear');
        
        % CONCAT AND SORT
        data_video_spk=sortrows([[SpikeFile X Y A VEL VELidx ones(size(SpikeFile,1),1)];[data_video_nospk,zeros(length(data_video_nospk),1)]],1);
        
        % VELO FILTER BASED ON INDEX CREATED ABOVE
        data_video_nospk(logical(in),:)=[];
        data_video_nospk(:,6)=[];
        nspkbefore=sum(data_video_spk(:,7)==1);
        data_video_spk(data_video_spk(:,6)==1,:)=[];
        nspkafter=sum(data_video_spk(:,7)==1);
        data_video_spk(:,6)=[];
        data_video_nospk(isnan(data_video_nospk(:,1)),:)=[];
        data_video_spk(isnan(data_video_spk(:,1)),:)=[];
        disp(['Velocity Filter: ',num2str(round(sum(in)/30)),' Seconds ', 'and ',num2str(nspkbefore-nspkafter),' out of ',num2str(length(SpikeFile)),' spikes'])
        
        % IFR
%         %ifr=IFR(data_video_spk(data_video_spk(:,6)==1,1),data_video_spk(:,6),round(length(data_video_nospk)/30),30);
%         ifr=IFR(data_video_spk,30);
%         VIFR=[data_video_spk(:,5),ifr];
%         VIFR(isnan(VIFR),:)=[];
%         SpeedScore=corr(VIFR(:,1),VIFR(:,2));
%         clear ifr VIFR
        SpeedScore=NaN;
%         % THETA MODULATION
%         [thetaindex,thetapeak,cor,~]=thetamodulation(SpikeFile./10000000);
%         data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).autocor{i,event}=cor;
       thetaindex=NaN;
        % OVERALL FIRING RATE
        OverallFR=(sum(data_video_spk(:,6))/size(data_video_nospk,1))*sampleRate;
        
        % NUMBER OF SPIKES
        nSpikes=sum(data_video_spk(:,6));
        
        % 6 degree bins
        da=pi/30;
        angBins=da/2:da:2*pi-da/2;
        % Occupancy
        histAng=hist(data_video_nospk(:,4),angBins);
        %creating variable to check if all directional bins were sampled. 
        if length(find(histAng))==60
            bins_Sampled=1;
        else
            bins_Sampled=0;
        end
        % Number of spikes per bin
        spkPerAng=hist(data_video_spk(data_video_spk(:,6)==1,4),angBins);
        
        % Tuning
        hdTuning=(spkPerAng./histAng)*sampleRate;
        % remove nan & inf
        hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
         data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).hdTuning{i,event}=hdTuning;
        
         %HD STATS =======================================================
         % calculates mean vector length, Directional information and half-width content for determination of
        %directionally modulated cells as well as basic characteristics. 
        mean_vector_length = circ_r(angBins',hdTuning',deg2rad(6));
        [peak_Firing_Rate,peakIdx]=max(hdTuning);
        preferred_Direction=peakIdx*6;
        
        [halfWidth,~]=FindHalfWidth(hdTuning);
        
        DIC= computeDIC(histAng,hdTuning,OverallFR);
        
        [within_Coeff,within] = within_HDstability(data_video_spk,data_video_nospk,sampleRate,4,6);
        data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).InterTrialStability{i,event}=within;

        %Coherence 
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
        
        clear binSize l hdTmp Npoints gw 
        
        %Preferred Direction drift
        [PDD_max,PDD_min, meanHD,N,plotIdx] = HDdrift(ExtractedAngle,data_video_spk,data_video_nospk,0);
        data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).drift.meanHD{i,event}=meanHD;
        data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).drift.binnedFR{i,event}=N;
        data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).drift.plotIdx{i,event}=plotIdx;

        
        %Create structured cell array so that shuffling can occur at
        %correct frame rate (accounts for interpolated timestamps that
        %alters frame rate of data aka data_video_spk. 
        secIdx=data_video_nospk(30:30:length(data_video_nospk),1);
        [I,row]=ismember(data_video_spk(:,1),secIdx);
        I=find(I);
        I=[1;I;length(data_video_spk)];
        
        for s=1:length(I)
            if s==length(I)
                datacell{s,1}=data_video_spk(I(s):I(end),:);
                break
            end
            datacell{s,1}=data_video_spk(I(s):I(s+1)-1,:);
        end
        

        %Shuffled R-Length (see methods Langston et al.,2010 & Boccara et
        %al., 2010) will collect shuffled data for cell and store in data
        %structure. Rlength for all shuffled data will then be assessed in
        %terms of populuation. 
        tempframes=data_video_spk;
        shuff_max=(sessionduration_sec-20);
        bispk=[];
        for ishuff=1:400
            shift=randi([20 round(shuff_max)]);
            tempcell=circshift(datacell,shift);
            for open=1:length(tempcell)
                bispk=[bispk;tempcell{open}(:,6)];
            end
            tempframes(:,6)=bispk;
            bispk=[];
            spkPerAng=hist(tempframes(tempframes(:,6)==1,4),angBins);
            hdTuning_shuff=(spkPerAng./histAng)*sampleRate;
            rlength_shuff(ishuff,1) = circ_r(angBins',hdTuning_shuff',deg2rad(6));
        end
         data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).rlength_shuff{i,event}=rlength_shuff;
        
        
        if  mean_vector_length >= prctile(rlength_shuff,95)
            pass_shuff=1;
        else
            pass_shuff=0;
        end
        
        %OTHER SPATIAL CHARACTERISTICS ====================================
        [SmoothRateMap,nBinsx,nBinsy,occ,Coherence]=bindata(data_video_nospk,sampleRate,data_video_spk(data_video_spk(:,6)==1,:),linear_track,track_length);
        data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ratemap{i,event}=SmoothRateMap;
        SmoothRateMap2=SmoothRateMap;
        [field,FieldWidth]=FindFF2D(SmoothRateMap2);
        FieldWidth=FieldWidth*(track_length/length(field));
        if isnan(field)==0
            nanPlacement=isnan(SmoothRateMap2);
            SmoothRateMap2(~field)=0;
            SmoothRateMap2(nanPlacement)=NaN;
            Map4bordertemp=SmoothRateMap2;
            Map4bordertemp(~field)=NaN;
            infieldFR=nanmean(reshape(Map4bordertemp,[],1));
            Map4bordertemp=SmoothRateMap;
            Map4bordertemp(logical(field))=NaN;
            outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
            clear Map4bordertemp
        else
            infieldFR=0;outfieldFR=0;
        end
        
       [~, filename] = fileparts(ID{i}); 
       tetrode=regexp(filename,'.','match');
       eegfile=cellstr(strcat(path,filesep,'CSC',tetrode(end),'.ncs'));
            
       [phaselockR,phaselockPval,ThetaStats]=EEGWorking_HD(eegfile{1},data_video_spk(data_video_spk(:,6)==1,:),StartofRec,EndofRec,event,track_length);
       clear tetrode eegfile normalizedD
        
        % CALCULATE IC FROM NSMA 
        rY=reshape(SmoothRateMap,nBinsx*nBinsy,1);rY(isnan(rY))=0;rY(isinf(rY))=0;
        occRSHP=reshape(occ,nBinsx*nBinsy,1);occSUM=sum(occRSHP);pX=occRSHP./occSUM;
        [nBins,nCells]=size(rY);relR=rY./kron(ones(nBins,1),pX'*rY);
        log_relR=log2(relR);log_relR(isinf(log_relR))=0;
        InformationContent=sum(kron(pX,ones(1,nCells)).*relR.*log_relR);
        
        %BORDER SCORE 
        [borderScore,Field2Wall]=BorderScore(SmoothRateMap2);
        Field2Wall=Field2Wall*(track_length/length(SmoothRateMap)); %F2wall converted to cm
        
        %PEAK RATE FOR PLACE FIELD 
        r_max=max(rY);
        PeakRate=r_max;

        % calculate sparsity from NSMA toolbox
        R=rY';if ~isa(R,'cell');R={R};end;RC=cat(2,R{:});[~,nCells]=size(RC);
        sparsity=sum(RC,2).^2./(nCells*sum(RC.^2,2));
        
        % calculate percent of active bins (estimate of field size - Royer et al., 2010)
        NumbActiveBins=numel(find(rY > r_max*0.20));
        clear R rY r_max RC log_relR nBins relR occRSHP occSUM pX nCells x_max
%         
%         measures=[bins_Sampled,mean_vector_length,halfWidth,DIC,peak_Firing_Rate,...
%             OverallFR,nSpikes,preferred_Direction,clusterquality(i,:),prop(i,:),...
%             pass_shuff,within_Coeff,Coherence,FieldWidth,infieldFR,outfieldFR,...
%             InformationContent,borderScore,Field2Wall,PeakRate,sparsity,NumbActiveBins,...
%             thetaindex,thetapeak,ThetaStats.MeanThetaFreq,ThetaStats.MeanOverallPow,ThetaStats.ThetaRatio,...
%             phaselockR,phaselockPval, HD_Coherence,PDD_max,PDD_min];

        measures=[bins_Sampled,mean_vector_length,halfWidth,DIC,peak_Firing_Rate,...
            OverallFR,nSpikes,preferred_Direction,clusterquality(i,:),prop(i,:),...
            pass_shuff,within_Coeff,Coherence,FieldWidth,infieldFR,outfieldFR,...
            InformationContent,borderScore,Field2Wall,PeakRate,sparsity,NumbActiveBins,...
            HD_Coherence,PDD_max,PDD_min];
        
        %SAVE SESSION DATA
        data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,event)=measures;
        
        clearvars -except  ExtractedAngle path data_video SpikeFile mclustpath tfile track_length...
            linear_track SmoothRateMap_Right SmoothRateMap_Left figures S...
            sampleRate i Overall_DirectionalityIndex StartofRec EndofRec event...
            timestamps vel_cmPerSec pixelDist grade data ratID video clusterquality prop ID sessionduration_sec
    end
end

% SAVE RESULTS TO MAT FILE
if exist('F:\P30_Recordings\data.mat','file')>0
    datatemp=data;
    load('F:\P30_Recordings\data.mat','data')
    data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')])=datatemp.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]);
end
save('D:\P30_Recordings\data.mat','data','-v7.3')
 wer
disp 'DONE'
