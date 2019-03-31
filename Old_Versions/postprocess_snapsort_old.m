
%%%%%%%% DATA EXTRACTION & POSTPROCESSING from Neuralynx SNAP sort file %%%%%%%%%%%%%
%
% THIS SCRIPT CREATES FIGURES AND STATS FOR PLACE & HD CELL DATA
% CAN HANDLE CIRCULAR ARENA & LINEAR TRACK DATA
%
% INPUTS: PATH TO RAW NON-PROCESSED DATA (NON-P FOLDER) / TRACK LENGTH / MAZE-TYPE / FIGURES
%
% OUTPUTS: RATE MAP / SPIKE ON PATH / TUNING CURVE / POLAR PLOT / R VS. L SPIKES / STATS / LFP ANALYSIS
% Ryan Harvey: 12/21/17

%    __    __    __    __    __    __    __    __    __    __    __
% __/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__
%   \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/


function postprocess_snapsort(path,track_length,linear_track,figures)
% close all;clc
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
clear com basedir

% sets path to snap folder & looks within for mat files
cd([path,filesep,'SNAPSorterResults']);
file=struct2table(dir( '**/*.mat'));
t=table2cell(file(:,1));
file=t(~contains(t,'_info.mat'),1);
info=t(contains(t,'_info.mat'),1);

% EXTRACT SPIKE TIMES, GRADES, CLUSTER QUALITY, ID, AND AVERAGE WAVEFORMS
% FROM SNAPSORT .MAT FILES
ns=1;grade=[];clusterquality=[];
for m=1:length(file)
    load(file{m});load(info{m})
    grade=[grade;final_grades'];
    clusterquality=[clusterquality;grades(:,1:10)];
    for n=1:length(unique(output(:,2)))
        ID{ns,1}=orig_filename;
        avgwave{ns,1}=means{1,n};
        S{ns,1}=output(output(:,2)==n,1); ns=ns+1;
    end
end
clear file info t file grades final_grades orig_filename output m confidence means ns n

disp(path)
if ~exist('S','var')
    disp('No Clusters... Skipping')
    return
end
disp(['    ',num2str(length(S)),' Cells'])

% CALCULATE WAVEFORM PROPERTIES
prop=waveformprop(avgwave);

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
varnames={'InformationContent','Coherence','Sparsity','PeakRate',...
    'OverallFiringRate','Field2Wall','FieldWidth','nSpikes','mean_vector_length',...
    'preferred_Direction','Direct_infoContent','DirectionalityIndex',...
    'PhPrecessSlope','PH_RSquared','lapPhPrecessSlope','lapPH_RSquared',...
    'PhlapCorrelation','lapPhaseOffset','PhcircLinCorr','Phpval','slopeCpU',...
    'phaseOffset','MeanThetaFreq','MeanThetaPow','Displacement',...
    'infieldFR','outfieldFR','Cluster Grade','borderScore','E','NumbActiveBins',...
    'DisplacementCorr','SpeedScore','IntraTrialStability','nspikes infield',...
    'ThetaPhaseLock R','ThetaPhaselock Pval','thetaindex','thetapeak',...
    'LRatio','Tightness','ShortISI','Incompleteness','IsolationDistance',...
    'NumberOfSpikes','StationInTime','TempMatch','BDistanceClust','BDistanceSpike'...
    'WaveformPeak','WaveformValley','spikewidth','halfwidth','peak2valleyratio','peak2valleyslope','SpikeAmplitude',...
    'temporalstability','nlaps','rateoverlap','fieldoverlap',...
    'lap_perm_stability','stabilityoverlaps','meanstability','spatialcorrelation'};

ratID=strsplit(path,filesep);
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).varnames=varnames;
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).Spikes=cellfun(@(x) x*1000000,S,'un',0);
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).avgwave=avgwave;

% ADD ENTIRE VIDEO TO DATA
[table,SumofNaN1,SumofNaN2,~]=NonDetects(video);

% if SumofNaN1>SumofNaN2; VelIndex=[4 5]; else; VelIndex=[2 3]; end
% [vel_cmPerSec,~,~] = InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,30);
% table = table((vel_cmPerSec<=100),:);

% fix tracker jumps
[table(:,2),table(:,3)]=trackerjumps(table(:,2),table(:,3));
[table(:,4),table(:,5)]=trackerjumps(table(:,4),table(:,5));

padDsData=[repmat(table(1,:),30,1); table; repmat(table(end,:),30,1)];
padDsData=[padDsData(:,1),smooth(padDsData(:,2),15,'lowess'),smooth(padDsData(:,3),15,'lowess'),smooth(padDsData(:,4),15,'lowess'),smooth(padDsData(:,5),15,'lowess')];
table=(padDsData(30+1:end-30,:));
data_video=double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2)]);
if isequal(linear_track,'yes')
    templinear=data_video(data_video(:,1)>StartofRec(1) & data_video(:,1)<EndofRec(1),:);
    [X,Y]=ParameterizeToLinearTrack2(templinear(:,2),templinear(:,3));
    data_video(data_video(:,1)>StartofRec(1) & data_video(:,1)<EndofRec(1),2:3)=[X Y];
end
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).frames=data_video;
clear data_video X Y padDsData vel_cmPerSec SumofNaN1 SumofNaN2 VelIndex templinear

% ADD EVENT TO DATA
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).events=[StartofRec;EndofRec];

% EXTRACT LFP FROM ALL TETRODES FOR LATER USE
% tetrodes=unique(extractBetween(ID,'TT','.ntt'));
% eegfile=strcat(repmat([path,filesep,'CSC'],length(tetrodes),1),tetrodes,repmat(['.ncs'],length(tetrodes),1));

channels=table2cell(struct2table(dir([path,filesep,'*.ncs'])));
eegfile=[strcat(channels(:,2),filesep,channels(:,1))];

[EEG_DownSampledTimestamps,EEG_DownSampledData,EEGthetaData]=downsampleLFP(eegfile);
clear tetrodes eegfile

% START OF EVENT LOOP 1:nSESSIONS
for event=1:length(StartofRec)
    if event==1 && isequal(linear_track,'yes');linear_track='yes';else;linear_track = 'no';end
    
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
    disp(['SESSION_',num2str(event),' WAS ',num2str(sessionduration),' MIN'])
    if sessionduration<2;disp('SESSION TOO SHORT...SKIPPING');continue;end
    clear sessionduration
    
    % DEAL WITH NON DETECTS
    [table,SumofNaN1,SumofNaN2,~]=NonDetects(table);
    
    % VELOCITY FILTER & SMOOTHING
    % Arena diameter
    ratID=strsplit(path,filesep);
    if isequal(linear_track,'no') && event==4
        track_length = 100;
    elseif isequal(linear_track,'no')
        track_length = 76.5;
    end
    
    %     CHOOSE LED WITH BEST TRACKING
    if SumofNaN1>SumofNaN2; VelIndex=[4 5]; else; VelIndex=[2 3]; end
    clear SumofNaN1 SumofNaN2
    
    % CALCULATE VEL TO REMOVE JUMPS IN DATA
    % fix tracker jumps
    [table(:,2),table(:,3)]=trackerjumps(table(:,2),table(:,3));
    [table(:,4),table(:,5)]=trackerjumps(table(:,4),table(:,5));
    
    
    % SMOOTHING RAW XY DATA FROM EACH LED
    padDsData=[repmat(table(1,:),30,1); table; repmat(table(end,:),30,1)]; %pad ends with last data point
    padDsData=[padDsData(:,1),smooth(padDsData(:,2),15,'lowess'),smooth(padDsData(:,3),15,'lowess'),smooth(padDsData(:,4),15,'lowess'),smooth(padDsData(:,5),15,'lowess')];
    table=(padDsData(30+1:end-30,:)); %remove pad
    clear padDsData
    
    % GET VELOCITY
    [vel_cmPerSec,vel_abs,pixelDist]=InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,sampleRate);
    clear VelIndex
    % SMOOTH VELOCITY WITH A WIDTH OF .8 SECONDS
    vel_cmPerSec=smooth(vel_cmPerSec,24);
    vel_abs=smooth(vel_abs,24);
    
    % CONCAT TIMESTAMP/XY POSITION/ANGLE TOGETHER
    data_video=double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2)]);
    
    % STRAIGHTEN LINEAR TRACK
    if isequal(linear_track,'yes')
        [X,Y]=ParameterizeToLinearTrack2(data_video(:,2),data_video(:,3));
        data_video=[data_video(:,1) X Y];
    end
    clear X Y
    
    % FINDING ANGLE (GETS ANGLE FROM RAW XY FOR EACH LED)
    ExtractedAngle=XYangle(data_video(:,2),data_video(:,3));
    data_video=[data_video,[ExtractedAngle(1);ExtractedAngle],[vel_cmPerSec(1);vel_cmPerSec]]; % remove first point to fit with vel
    
    % EXTRACT BASIC MOVEMENT DATA
    BasicLoco.AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*sampleRate));
    BasicLoco.OverallDistTraveled=sum(vel_abs*pixelDist);
    BasicLoco.MeanVelocity=mean(vel_cmPerSec);
    data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).BasicLoco=BasicLoco;
    
    clear vel_cmPerSec vel_abs ExtractedAngle BasicLoco
    
    % START OF MAIN SPIKE TO PATH AND DATA ANALYSIS LOOP
    for i=1:length(S)
        cell=strsplit(ID{i},filesep);
        disp([cell{end},'   Cell: ',num2str(i)])
        
        SpikeFile=S{i}*1000000;
        
        % RESTRICT DATA BY START AND END OF EVENT
        if exist('timestamps','var')==1
            SpikeFile=SpikeFile(SpikeFile(:,1)>StartofRec(event) & SpikeFile(:,1)<EndofRec(event),:);
        end
        
        % FIND INDEX FOR VELOCITY FILTER (<2cm/second for 2 second)
        
        in=contiguousframes(data_video(:,5)<2,60);
        data_video_nospk=[data_video,in];
        
        % INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
        X=interp1(data_video_nospk(:,1),data_video_nospk(:,2),SpikeFile,'linear');
        Y=interp1(data_video_nospk(:,1),data_video_nospk(:,3),SpikeFile,'linear');
        A=interp1(data_video_nospk(:,1),data_video_nospk(:,4),SpikeFile,'linear');
        VEL=interp1(data_video_nospk(:,1),data_video_nospk(:,5),SpikeFile,'linear');
        VELidx=interp1(data_video_nospk(:,1),data_video_nospk(:,6),SpikeFile,'linear');
        
        % CONCAT AND SORT
        data_video_spk=sortrows([[SpikeFile X Y A VEL VELidx ones(size(SpikeFile,1),1)];[data_video_nospk,zeros(length(data_video_nospk),1)]],1);
        clear TS X Y A VEL VELidx
        
        % VELO FILTER BASED ON INDEX CREATED ABOVE
        data_video_nospk(logical(in),:)=[];
        data_video_nospk(:,6)=[];
        %         nspkbefore=sum(data_video_spk(:,7)==1);
        data_video_spk(data_video_spk(:,6)==1,:)=[];
        %         nspkafter=sum(data_video_spk(:,7)==1);
        data_video_spk(:,6)=[];
        
        %         disp(['Velocity Filter: ',num2str(round(sum(in)/30)),' Seconds ', 'and ',num2str(nspkbefore-nspkafter),' out of ',num2str(length(SpikeFile)),' spikes'])
        clear in nspkbefore nspkafter
        
        % SEE IF SPIKES REGULARLY OCCUR OVER TIME
        temporal_stability=temporalstability(data_video_spk(data_video_spk(:,6)==1,1),data_video_spk(:,1));
        
        % calculate IFR
        ifr=IFR(data_video_spk,30);
        % calculate running speed over 200ms bins
        n = 6; % average every n values
        b = arrayfun(@(i) mean(data_video_nospk(i:i+n-1,5)),1:n:length(data_video_nospk)-n+1)'; % the averaged vector
        % correlation ifr to average running speed to see if firing rate
        % increases as a function of running speed
        
        % make sure two vectors are same length
        if length(ifr)~=length(b)
            [~,I]=max([length(ifr);length(b)]);
            if I==1
                ifr(length(b)+1:end)=[];
            elseif I==2
                b(length(ifr)+1:end)=[];
            end
        end
        SpeedScore=corr(ifr,b);
        clear ifr b n I
        
        % SPLIT UP DATA BY RIGHT AND LEFT RUNS FOR LINEAR TRACK DATA
        if isequal(linear_track,'yes')
            runiteration=2;
            [right,left,DirectionalityIndex,Displacement,nlaps,rateoverlap,fieldoverlap,spatialcorrelation]=RightVsLeft(data_video_nospk,data_video_spk,track_length,sampleRate);
        else
            runiteration=1;
        end
        data_video_smoothfilt3=data_video_spk;
        for iruns=1:runiteration
            if isequal(linear_track, 'yes') % UNPACK STRUCTURE
                if iruns==1
                    occ=right.occ; nBinsx=right.nBinsx; nBinsy=right.nBinsy;
                    Coherence=right.Coherence;
                    SmoothRateMap=right.SmoothRateMap;
                    data_video_smoothfilttemp=right.dataspks;
                    lap_perm_stability=right.lap_perm_stabilityR;
                    stabilityoverlaps=right.stabilityoverlapsR;
                    meanstability=right.meanstabilityR;
                elseif iruns==2
                    occ=left.occ; nBinsx=left.nBinsx; nBinsy=left.nBinsy;
                    Coherence=left.Coherence;
                    SmoothRateMap=left.SmoothRateMap;
                    data_video_smoothfilttemp=left.dataspks;
                    lap_perm_stability=left.lap_perm_stabilityL;
                    stabilityoverlaps=left.stabilityoverlapsL;
                    meanstability=left.meanstabilityL;                   
                end
                
                [~,ia]=setdiff(data_video_smoothfilt3(:,1),data_video_smoothfilttemp(:,1));
                plotframes=data_video_smoothfilt3;
                plotframes(ia,:)=NaN;
                
                data_video_spk=data_video_smoothfilttemp;
                clear ia data_video_smoothfilttemp
            end
            spks_VEL=data_video_spk(data_video_spk(:,6)==1,:);
            
            % THETA MODULATION
            % convert ts to secs
            sec=linspace(0,(length(data_video_nospk(:,1))/sampleRate),length(data_video_nospk(:,1)));
            spktssec=interp1(data_video_nospk(:,1),sec',spks_VEL(:,1));
            % calc theta mod with
            [thetaindex,thetapeak,cor,~]=thetamodulation(spktssec);
            clear spktssec
            
            %%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%
            % BIN & SMOOTH CIRCLE DATA
            if isequal(linear_track,'no')
                [SmoothRateMap,nBinsx,nBinsy,occ,Coherence]=bindata(data_video_nospk,sampleRate,spks_VEL,linear_track,track_length);
                SmoothRateMap2=SmoothRateMap;
                [field,FieldWidth]=FindFF2D(SmoothRateMap2);
                FieldWidth=FieldWidth*(track_length/length(field));
                if isnan(field)==0
                    % in field / out field FR
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
                    % spikes in field
                    [row,col]=find(field);k=boundary(row,col);bound=[col(k),row(k)];
                    if length(bound)>2
                        frames=[rescale(data_video_spk(:,2),1,length(field)),rescale(data_video_spk(:,3),1,length(field))];
                        in=inpolygon(frames(:,1),frames(:,2),bound(:,1),bound(:,2));
                        nspikes_infield=sum(data_video_spk(in & data_video_spk(:,6)==1,6));
                        clear frames in 
                    else
                        nspikes_infield=0;
                    end
                    clear row col k bound
                else
                    infieldFR=0;outfieldFR=0;nspikes_infield=0;
                end
                
                % CALCULATE LANDMARK CONTROL
                if event==3
                    [Displacement,DisplacementCorr]=Displacement2(SmoothRateMap,data.(ratID{end-1}).(['S',...
                        strjoin(regexp(ratID{end},'\d*','Match'),'')]).ratemap.(['Cell',num2str(i)]).(['session',num2str(event)]));
                end
                
                [borderScore,Field2Wall]=BorderScore(SmoothRateMap2);
                Field2Wall=Field2Wall*(track_length/length(SmoothRateMap));
                
                % find field
                [row,col]=find(field);k=boundary(row,col);bound=[col(k),row(k)];
                if length(bound)>2
                    
                    % Eccentricity
                    E=Eccentricity(col(k),row(k));
                    
                    % SET UP PASSES THROUGH FIELD FOR PHASE PRECESSION
                    % convert ts to ms
                    ms=[linspace(0,(length(data_video_nospk)/sampleRate),length(data_video_nospk))]*1000;
                    ms=interp1(data_video_nospk(:,1),ms,data_video_spk(:,1));
                    % rescale xy to size of rate map so you can find xy cors in field
                    frames=[rescale(data_video_spk(:,2),1,length(SmoothRateMap)),rescale(data_video_spk(:,3),1,length(SmoothRateMap))];
                    % find in durations field
                    in=inpolygon(frames(:,1),frames(:,2),bound(:,1),bound(:,2));
                    % pull out each pass and store in cell array
                    dsig=diff([0 (abs(in')>=eps) 0]);
                    startIndex=find(dsig>0);endIndex=find(dsig<0)-1;
                    stringIndex=(endIndex-startIndex+1>=0);
                    startIndex=startIndex(stringIndex);
                    endIndex=endIndex(stringIndex);
                    indices=zeros(1,max(endIndex)+1);
                    indices(startIndex)=1;
                    indices(endIndex+1)=indices(endIndex+1)-1;
                    indices=find(cumsum(indices));
                    in=zeros(length(in),1);in(indices',1)=1;
                    pass=find(diff([0 indices])>1==1);
                    if isempty(pass) || length(pass)==1
                        spks_VEL4LFP=NaN;
                        occ4Ph=NaN;
                        fieldbound=NaN;
                    else
                        if pass(1)~=1
                            pass=[1,pass];
                        end
                        for j=1:length(pass)
                            if j+1>length(pass);break;end
                            tempidx=indices(pass(j):pass(j+1)-1);
                            storets{j,1}=[ms(tempidx),data_video_spk(tempidx,1),...
                                frames(tempidx,:),data_video_spk(tempidx,4:end)];
                        end
                        % CRITERIA FOR PASSES THROUGH FIELD:
                        %   FILTER OUT < average 3CM/SEC per pass
                        %   FILTER OUT LESS THAN 200MS
                        %   FILTER OUT LESS THAN 5 SPIKES OVER AT LEAST 4 THETA CYCLES
                        %   MAX ISI <1000ms
                        x=[];
                        for j=1:length(storets)
                            temp=storets{j,1};
                            if mean(temp(:,6))>3 && (temp(end,1)-temp(1,1))>200 ...
                                    && sum(temp(:,end))>20 && max(diff(temp(:,1)))<1000
                                % make x linear
                                x=[x;linspace(0,1,size(temp,1))'];
                                continue
                            else
                                storets{j,1}=[];
                            end
                        end
                        
                        % combine cell array of passes and create spike matrix
                        % & occ matrix for use in PHprecession
                        tempframes=vertcat(storets{:});
                        if isempty(tempframes)
                            spks_VEL4LFP=NaN;
                            occ4Ph=NaN;
                            fieldbound=NaN;
                        else
                            spks_VEL4LFP=tempframes(:,1)/1000;
                            tempframes(:,3)=x;
                            occ4Ph=tempframes(tempframes(:,end)==0,[1,3:end]);
                            occ4Ph(:,1)=occ4Ph(:,1)/1000;
                            occ4Ph(:,3)=zeros(length(occ4Ph),1);
                            fieldbound=[0 1];
                        end
                        clear x tempframes storets temp j ms frames tempidx in dsig indices startIndex endIndex
                    end
                else
                    spks_VEL4LFP=data_video_spk(data_video_spk(:,6)==1,:);
                    occ4Ph=data_video_spk(data_video_spk(:,6)==0,:);
                    E=NaN;
                    fieldbound=NaN;
                end
            else
                % LINEAR
                [~,FieldWidth,field]=FindFF(SmoothRateMap,data_video_spk,pixelDist,track_length);
                spks_VEL4LFP=data_video_spk(data_video_spk(:,5)>5 & data_video_spk(:,6)==1,:);
                
                % base field boundaries on velocity filted data
                [~,~,fieldvelofilt]=FindFF(SmoothRateMap,data_video_spk(data_video_spk(:,5)>5,:),pixelDist,track_length);
                fieldidx=find(fieldvelofilt);
                temp=rescale([1:length(fieldvelofilt), min(fieldidx),max(fieldidx)],0,1);
                fieldbound=temp(end-1:end);
                clear fieldvelofilt temp
                
                ts=sort(data_video_nospk(:,1));
                sec=linspace(0,(length(ts)/sampleRate),length(ts));
                spks_VEL4LFP(:,1)=interp1(ts',sec',spks_VEL4LFP(:,1));
                
                occ4Ph=[interp1(ts,sec,data_video_spk(:,1)),data_video_spk(:,2:end)];
                occ4Ph=occ4Ph(occ4Ph(:,6)==0,:);
                
                % in field / out field FR
                Map4bordertemp=SmoothRateMap;Map4bordertemp(~field)=NaN;
                infieldFR=nanmean(reshape(Map4bordertemp,[],1));if isnan(infieldFR)==1;infieldFR=0;end
                Map4bordertemp=SmoothRateMap;Map4bordertemp(logical(field))=NaN;
                outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
                clear Map4bordertemp
                
                % spikes in field
                [nspikes_infield,~,~]=FindFF(SmoothRateMap,data_video_spk,pixelDist,track_length);
                nspikes_infield=sum(nspikes_infield(:,end));
            end
            
            % LFP ANALYSIS
            % Create name of CSC file from cell and path name
            [~, filename] = fileparts(ID{i}); tetrode=regexp(filename,'.','match');
            % limit lfp by event
            trodeID=str2double(tetrode(end));
            if (length(StartofRec) + length(EndofRec))>2 
                samples=EEG_DownSampledData(trodeID,EEG_DownSampledTimestamps>StartofRec(event) & EEG_DownSampledTimestamps<EndofRec(event));
                ts=EEG_DownSampledTimestamps(EEG_DownSampledTimestamps>StartofRec(event) & EEG_DownSampledTimestamps<EndofRec(event));
                thetaphase=EEGthetaData(trodeID,EEG_DownSampledTimestamps>StartofRec(event) & EEG_DownSampledTimestamps<EndofRec(event));
            else
                samples=EEG_DownSampledData(trodeID,:);
                ts=EEG_DownSampledTimestamps;
                thetaphase=EEGthetaData(trodeID,:);
            end
            % calc phase precession stats
            [ThPrecess]=PHprecession(thetaphase,ts,spks_VEL4LFP,occ4Ph,fieldbound);
            % theta power
            [MeanThetaFreq,MeanOverallPow]=meanfreq(thetaphase,1000);
            
            clear tetrode eegfile ts samples thetaphase
            
            % calculate information content from nsma
            rY=reshape(SmoothRateMap,nBinsx*nBinsy,1);rY(isnan(rY))=0;rY(isinf(rY))=0;
            occRSHP=reshape(occ,nBinsx*nBinsy,1);occSUM=sum(occRSHP);pX=occRSHP./occSUM;
            [nBins,nCells]=size(rY);relR=rY./kron(ones(nBins,1),pX'*rY);
            log_relR=log2(relR);log_relR(isinf(log_relR))=0;
            InformationContent=sum(kron(pX,ones(1,nCells)).*relR.*log_relR);
            
            % Calculate the distance from place field to wall - peaks method
            r_max=max(rY);
            PeakRate=r_max;
            if isequal(linear_track,'yes')
                [~, x_max]=find(SmoothRateMap==r_max);
                Field2Wall=min([nBinsx-x_max x_max-0])*(track_length/length(SmoothRateMap));
            end
            
            % OVERALL FIRING RATE
            OverallFR=(size(spks_VEL,1)/size(data_video_nospk,1))*sampleRate;
            
            % NUMBER OF SPIKES
            nSpikes=size(spks_VEL,1);
            
            % calculate sparsity from NSMA toolbox
            R=rY';if ~isa(R,'cell');R={R};end;RC=cat(2,R{:});[~,nCells]=size(RC);
            sparsity=sum(RC,2).^2./(nCells*sum(RC.^2,2));
            
            % calculate percent of active bins (estimate of field size - Royer et al., 2010)
            NumbActiveBins=numel(find(rY > r_max*0.20));
            clear R rY r_max RC log_relR nBins relR occRSHP occSUM pX nCells x_max
            
            % GRID CELL ANALYSIS
            % GRID SCORE
            %             if isequal(linear_track,'no');
            %                 autocorr=SpatialAutoCorr(SmoothRateMap,length(SmoothRateMap));
            %                 gridout=GridScore_Sinusoidal(autocorr,length(SmoothRateMap));
            %                 gridscore=gridout.maxSinuGrid;
            %             end
            
            % INTRA-TRIAL STABILITY
            IntraTrialR=IntraTrialStability(data_video_spk,linear_track,track_length);
            
            % SPIKE DIRECTION (from Shawn Winter's code)
            [mean_vector_length,peak_Firing_Rate,preferred_Direction,~,...
                ~,Direct_infoContent,BinsNbSpikes,~,BinsAngle]...
                = HDCell(data_video_spk(:,6),data_video_spk(:,4),sampleRate);
            
            % -------------------------------CREATING PLOTS----------------------------
            if figures==1
                if isequal(linear_track, 'yes')
                    % plot spike on path
                    %                     fig=figure;subplot(5,1,1),plot(plotframes(:,2),plotframes(:,3),'LineWidth',1,'color','k');
                    %                     hold on;box off; axis off
                    %                     scatter(spks_VEL(:,2),spks_VEL(:,3),20,'filled','r');
                    %                     title([cell{end},' Cell: ',num2str(i),'  nSpikes: ',num2str(nSpikes)]);
                    %
                    fig=figure;subplot(4,1,1),plot(plotframes(:,2),plotframes(:,1),'LineWidth',1,'color','k');
                    hold on;box off; axis off
                    scatter(spks_VEL(:,2),spks_VEL(:,1),20,'filled','r');
                    title([cell{end},' Cell: ',num2str(i),'  nSpikes: ',num2str(nSpikes)]);
                    
                    % plot filled and smoothed rate map
                    imAlpha=ones(size([SmoothRateMap;SmoothRateMap]));
                    imAlpha(isnan([SmoothRateMap;SmoothRateMap]))=0;
                    subplot(4,1,2); imagesc([SmoothRateMap;SmoothRateMap],'AlphaData',imAlpha);
                    hold on; shading interp; colormap jet; axis off; box off;axis image
                    title(['InfoContent: ',num2str(round(InformationContent,3)),' F2W: ',num2str(round(Field2Wall,3))]);
                    
                    subplot(4,1,3), area(SmoothRateMap(1,:),'LineWidth',2,'EdgeColor',[0,0,0]+0.4,'FaceColor',[0,0,0]+0.8);
                    box off; xlim([1 nBinsx]); hold on
                    % plot field location
                    if sum(field)>0
                        FL=find(field==1);
                        plot([FL(1);FL(1)],[PeakRate;0],'LineWidth', 2, 'color','r')
                        plot([FL(end);FL(end)],[PeakRate;0],'LineWidth', 2, 'color','r')
                    end
                    title(['PR: ',num2str(round(PeakRate,3)),' OFR: ',num2str(round(OverallFR,3))]);
                    set(fig,'Position',[842 42 838 924]);
                    
                    % PLOT SCATTERED PHASE PRECESSION
                    if ~isnan(ThPrecess.scatteredPH)
                        subplot(4,1,4)
                        % rescale x distance for plotting
                        %                         fieldidx=find(field_restrict);
                        if ~isempty(fieldidx)
                            x=rescale(ThPrecess.scatteredPH(:,1),min(fieldidx),max(fieldidx));
                        else
                            x=NaN;
                        end
                        plot(x,ThPrecess.scatteredPH(:,2),'k.');hold on;
                        plot(x,ThPrecess.scatteredPH(:,2)+360,'r.')
                        ylim([0 720]); xlim([1 60])
                        set(gca, 'YTick', [0;240;480;720],'Box','off');
                        title(['Precess R^2: ',num2str(round(ThPrecess.RSquared,3)),', Slope: ',...
                            num2str(round(ThPrecess.slope,3))])
                        
                        % PLOT SMOOTHED RATEMAP
                        %                         subplot(5,1,5);
                        %                         h=pcolor(ThPrecess.smoothedPHmap);
                        %                         shading interp; colormap jet; hold on; box off; axis off; set(h, 'EdgeColor', 'none');
                        %                         title(['DOM: ',num2str(round(ThPrecess.DOM,3))])
                        
                        % PLOT AUTOCORR
                        %                         plot([linspace(-500,-1,100),0,linspace(1,500,100)], cor,'k'); hold on
                        %
                        %                         set(gca,'YTickLabel',[],'YTick',[],'XMinorTick','on','YMinorTick','off','LineWidth',1)
                        %                         line([0 0], ylim, 'linestyle', ':', 'color', [.7 .7 .7]);
                        % %                         set(gca, 'fontsize', 20, 'box', 'off');
                        %                         title(['Theta Index= ',num2str(thetaindex(j)),' Freq= ',num2str(peak)])
                    end
                    hold off
                    
                    % FOR CIRCULAR ARENA
                elseif isequal(linear_track,'no')
                    % plot spike on path
                    frames=[rescale(data_video_spk(:,2),1,length(SmoothRateMap)),rescale(data_video_spk(:,3),1,length(SmoothRateMap))];
                    spkfram=frames(data_video_spk(:,6)==1,:);
                    fig=figure;subplot(3,2,1);
                    plot(frames(:,1), frames(:,2),'LineWidth',1,'color','k');
                    hold on; axis off
                    scatter(spkfram(:,1), spkfram(:,2), 35, 'filled', 'r');
                    box off; axis image
                    plot(bound(:,1),bound(:,2),'LineWidth', 3, 'color', [.3 .3 .3])
                    title(['nSpikes: ',num2str(nSpikes),' F2W: ',num2str(round(Field2Wall,3))]);
                    
                    imAlpha=ones(size(SmoothRateMap));
                    imAlpha(isnan(SmoothRateMap))=0;
                    subplot(3,2,2);imagesc(SmoothRateMap,'AlphaData',imAlpha);
                    axis xy; colormap jet; axis off; hold on; box off; axis image;
                    title(['IC: ',num2str(round(InformationContent,3)),' PR: ',num2str(round(PeakRate,3)),' OFR: ',num2str(round(OverallFR,3))]);
                    
                    %smooth spike data
                    [pdf,~]=circ_ksdensity(spks_VEL(:,4), 0:359,'msni');
                    % Plot Smoothed firing rate x HEAD DIRECTION
                    subplot(3,2,3); plot(pdf,'LineWidth',2,'color','k')
                    axis tight; hold on; xlim ([0 360]); box off
                    title(['PR: ',num2str(round(peak_Firing_Rate,3)),' D IC: ',num2str(round(Direct_infoContent,3))])
                    
                    % Firing rate x HD polar plot for the nonsmoothed data above
                    subplot(3,2,4); polarplot = polar(BinsAngle([1:60 1]),BinsNbSpikes([1:60 1]),'b');
                    set(polarplot,'linewidth',3,'color','k'); axis off
                    title(['RLength: ',num2str(round(mean_vector_length,3)),' Pref Dir: ',num2str(round(preferred_Direction,3))]);
                    set(0,'Showhiddenhandles','on')
                    set(fig,'Position',[842 42 838 924]);
                    %                     fig1 = figure(Session);
                    
                    if ~isnan(ThPrecess.scatteredPH)
                        subplot(3,2,5); plot(ThPrecess.scatteredPH(:,1),ThPrecess.scatteredPH(:,2),'k.');hold on;
                        plot(ThPrecess.scatteredPH(:,1),ThPrecess.scatteredPH(:,2)+360,'r.')
                        ylim([0 720]); xlim([0 1])
                        set(gca, 'YTick', [0;240;480;720],'Box','off');
                        title(['Precess',' R^2: ',num2str(round(ThPrecess.RSquared,3)),', Slope: ',...
                            num2str(round(ThPrecess.slope,3)),', Corr: ',num2str(round(ThPrecess.Correlation,3))])
                        hold off
                        
                        % PLOT SMOOTHED RATEMAP
                        subplot(3,2,6); h=pcolor(ThPrecess.smoothedPHmap);
                        shading interp; colormap jet; hold on; box off; axis off; set(h, 'EdgeColor', 'none');
                        title(['DOM: ',num2str(round(ThPrecess.DOM,3))])
                    end
                    hold off
                end
                pause(.0000001)
            end
            close all
            
            % pack results for NaN measures given maze type
            if isequal(linear_track,'no')
                DirectionalityIndex=NaN;
                nlaps=NaN;
                rateoverlap=NaN;
                fieldoverlap=NaN;
                lap_perm_stability=NaN;
                stabilityoverlaps=NaN;
                meanstability=NaN;
                spatialcorrelation=NaN;
            else
                borderScore=NaN; E=NaN;DisplacementCorr=NaN;
            end
            if ~exist('Displacement','var');Displacement=NaN;DisplacementCorr=NaN;end
            
            measures=[InformationContent,Coherence,sparsity,PeakRate,...
                OverallFR,Field2Wall,FieldWidth,nSpikes,mean_vector_length,...
                preferred_Direction(1),Direct_infoContent,DirectionalityIndex,...
                ThPrecess.slope,ThPrecess.RSquared,ThPrecess.lapSlope,ThPrecess.lapR2,...
                ThPrecess.lapCorrelation,ThPrecess.lapPhaseOffset,ThPrecess.circLinCorr,...
                ThPrecess.pval,ThPrecess.slopeCpU,ThPrecess.phaseOffset,...
                MeanThetaFreq,MeanOverallPow,...
                Displacement,infieldFR,outfieldFR,grade(i),borderScore,E,NumbActiveBins,...
                DisplacementCorr,SpeedScore,IntraTrialR,nspikes_infield,...
                ThPrecess.phaselock.Rlength,ThPrecess.phaselock.Pval,thetaindex,thetapeak,...
                clusterquality(i,:),prop(i,:),temporal_stability,nlaps,rateoverlap,fieldoverlap,...
                lap_perm_stability,stabilityoverlaps,meanstability,spatialcorrelation];
                        
            if event==1 && iruns==1
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,event)=measures;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ratemap.(['Cell',num2str(i)]).(['session',num2str(event)])=SmoothRateMap;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).thetaautocorr.(['Cell',num2str(i)]).(['session',num2str(event)])=cor;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ThPrecess{i,event}=ThPrecess.scatteredPH;
                %                 data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).frames.(['Cell',num2str(i)]).(['session',num2str(event)])=data_video_smoothfilt;
            elseif event==1 && iruns>1
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,event+1)=measures;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ratemap.(['Cell',num2str(i)]).(['session',num2str(event+1)])=SmoothRateMap;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).thetaautocorr.(['Cell',num2str(i)]).(['session',num2str(event+1)])=cor;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ThPrecess{i,event+1}=ThPrecess.scatteredPH;
                
                %                 data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).frames.(['Cell',num2str(i)]).(['session',num2str(event+1)])=data_video_smoothfilt;
            elseif event>1
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,event+1)=measures;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ratemap.(['Cell',num2str(i)]).(['session',num2str(event+1)])=SmoothRateMap;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).thetaautocorr.(['Cell',num2str(i)]).(['session',num2str(event+1)])=cor;
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ThPrecess{i,event+1}=ThPrecess.scatteredPH;
                
                %                 data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).frames.(['Cell',num2str(i)]).(['session',num2str(event+1)])=data_video_smoothfilt;
            end
        end
        clearvars -except path data_video SpikeFile mclustpath tfile track_length...
            linear_track SmoothRateMap_Right SmoothRateMap_Left figures S...
            sampleRate i Overall_DirectionalityIndex StartofRec EndofRec event...
            timestamps vel_cmPerSec pixelDist grade data ratID video clusterquality prop ID...
            EEG_DownSampledTimestamps EEGthetaData EEG_DownSampledData
    end
end

% SAVE RESULTS TO MAT FILE
% if exist('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','file')>0
%     datatemp=data;
%     load('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','data')
%     data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')])=datatemp.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]);
% end
% save('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','data','-v7.3')

% save mat file to temp file
save(['D:\Place_Cell_Data\RawPAE_PlaceCell\temp\',(ratID{end-1}),'_',(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')])],'data','-v7.3')


disp 'DONE'

