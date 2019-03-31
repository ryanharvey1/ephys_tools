% High level function that processes single unit and lfp data
%
% output variables listed in variable 'varnames' below
%
% Ryan E. Harvey
%
function data=postprocess(path,track_length,linear_track,figures)
com=which('postprocess');
com=strsplit(com,filesep);

basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,'BClarkToolbox',filesep,'Analyses',filesep,'spikeCode',filesep,'MEX'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analyses',filesep,'spikeCode'],...
    [basedir,filesep,'chronux_2_11'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Utils'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Utils',filesep,'Inpaint_nans'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'Visualize'],...
    [basedir,filesep,'BClarkToolbox',filesep,'Analysis',filesep,'LFP'],...
    [basedir,filesep,'CircStat2012a'],...
    [basedir,filesep,'FMAToolbox'],...
    genpath([basedir,filesep,'neuralynxmex']),...
    [basedir,filesep,'BClarkToolbox',filesep,'Analyses',filesep,'spikeCode',filesep,'read_nvt'],...
    [basedir,filesep,'Pass_Index']);

clear com basedir

% sets path to snap folder & looks within for mat files
cd([path,filesep,'SNAPSorterResults']);
file=struct2table(dir( '**/*.mat'));
t=table2cell(file(:,1));
file=t(~contains(t,'_info.mat'),1);
info=strcat(erase(file,'.mat'),repmat('_info.mat',length(file),1));

% EXTRACT SPIKE TIMES, GRADES, CLUSTER QUALITY, ID, AND AVERAGE WAVEFORMS
% FROM SNAPSORT .MAT FILES
ns=1;grade=[];clusterquality=[];
for m=1:length(file)
    load(file{m});load(info{m})
    grade=[grade;final_grades'];
    clusterquality=[clusterquality;grades(:,1:10)];
    cells=unique(output(:,2));
    idxmeans=1;
    for n=1:length(cells)
        tetrode{ns,1}=file{m};
        cellnum(ns,1)=cells(n);
        ID{ns,1}=orig_filename;
        avgwave{ns,1}=means{1,idxmeans};
        S{ns,1}=output(output(:,2)==cells(n),1);
        ns=ns+1;
        idxmeans=idxmeans+1;
    end
end

disp(path)
if ~exist('S','var')
    disp('No Clusters... Skipping')
    return
end
disp(['    ',num2str(length(S)),' Cells'])

% Snap sorter outputs ts in ms. We need to convert to microseconds
for spikes=1:length(S)
    if sum(S{spikes}<1e8)==length(S{spikes})
        S{spikes}=S{spikes}*1000000;
    end
end
clear file info t file grades final_grades orig_filename output m confidence means ns n spikes

% CALCULATE WAVEFORM PROPERTIES
prop=waveformprop(avgwave);

% READ VIDEO DATA
[ts,x,y,angles] = Nlx2MatVT([path,filesep,'VT1.nvt'],[1,1,1,1,0,0],0,1);

% calculate sample rate from ts
tempts=(ts-ts(1))/10^6;
idx=find(tempts>=1);
sampleRate=idx(1)-1;
clear tempts idx
data.samplerate=sampleRate;

% DURATION OF SESSION (MIN)
if ((length(ts)/sampleRate)/60)<4
    disp('SESSION TOO SHORT...SKIPPING');
    return
end

ratID=strsplit(path,filesep);
data.rat=ratID{end-1};
data.sessionID=['S',strjoin(regexp(ratID{end},'\d*','Match'),'')];
data.date_processed=date;
data.session_path=path;

% SPLIT UP SESSIONS BY EVENT FILES
[~,StartofRec,EndofRec]=EventSplit(path);
if StartofRec(1)==1 && EndofRec(1)==1
    StartofRec=ts(1);
    EndofRec=ts(end);
end

% MAZE TYPES 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for mt=1:length(StartofRec)
    if mt<5 && contains(path,'ClarkP30')
        mazetype{mt}='Cylinder';
    elseif mt>=5 && contains(path,'ClarkP30')
        mazetype{mt}=input(['Maze Type for Session: ',num2str(mt),'  '],'s');
    end
    if contains(linear_track,'yes') && ~contains(path,'HPCatn') && mt==1
        mazetype{mt}='LinearTrack';
    elseif contains(linear_track,'yes') && contains(path,'HPCatn') && mt==1
        mazetype{mt}='CircularTrack';
    elseif contains(linear_track,'yes') && ~contains(path,'HPCatn') && mt>1
        mazetype{mt}='Cylinder';
    elseif contains(linear_track,'yes') && contains(path,'HPCatn') && mt>1
        mazetype{mt}='Box';
    end
end
data.mazetypes=mazetype; clear mt mazetypes

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
    'lap_perm_stability','stabilityoverlaps','meanstability','spatialcorrelation',...
    'egomod','bordermod','burstIdx'};

data.varnames=varnames;
data.Spikes=S;
data.spikesID.paths=ID;
data.spikesID.TetrodeNum=tetrode;
data.spikesID.CellNum=cellnum;
data.avgwave=avgwave;
clear varnames S avgwave

% CHECK FOR MANUAL COORDINATE CORRECTIONS
if exist([path,filesep,'restrictxy.mat'],'file')
    % run the following line if tracker errors from unplugs are present
%     manual_trackerjumps(ts,x,y,StartofRec,EndofRec,path);
    load([path,filesep,'restrictxy.mat'],'in')
    x(in==0)=NaN;
    y(in==0)=NaN;
    clear in
end

% FIX NON-DETECTS
[xtemp,ytemp]=FixPos(x',y',ts',round(0.1667*sampleRate));
tempangle=wrapTo360(fixNLXangle(angles',round(0.1667*sampleRate)))'-90;
tempangle(tempangle<0)=tempangle(tempangle<0)+360;

data_video=[ts',xtemp,ytemp,tempangle]; clear xtemp ytemp ts

if isequal(linear_track,'yes')
    templinear=data_video(data_video(:,1)>=StartofRec(1) & data_video(:,1)<=EndofRec(1),:);
    data.linear_track.nonlinearFrames=templinear;
    if sum(contains(ratID,'HPCatn'))>0
        [X,Y]=ParamToHorseshoeTrack(templinear(:,2),templinear(:,3));
    else
        [X,Y]=ParameterizeToLinearTrack2(templinear(:,2),templinear(:,3));
    end
    data_video(data_video(:,1)>=StartofRec(1) & data_video(:,1)<=EndofRec(1),2:3)=[X Y];    
end
if sum(contains(path,'HPCatn'))==1 && length(StartofRec)>1
    % restrict movement to inside the maze
    if exist([path,filesep,'data_video.mat'],'file')>0
        load([path,filesep,'data_video.mat'])
    end
end

% FIND MAX DIM OF BIGGEST MAZE TO CALCULATE VELOCITY 
if sum(contains(path,'HPCatn'))>0 
    track_length=360;
elseif sum(contains(path,'HPCatn'))==0 && size(StartofRec,2)==1 && track_length~=90
    track_length=120;
elseif contains(linear_track,'no')
    track_length=76.5;
end

% GET VELOCITY
[vel_cmPerSec,~,pixelDist]=InstaVel([data_video(:,2),data_video(:,3)],linear_track,track_length,sampleRate);

% vel_abs=smooth(vel_abs,24);
data_video=[data_video,[vel_cmPerSec(1);vel_cmPerSec]];

data.frames=data_video;
clear data_video X Y vel_cmPerSec templinear XY x y ts angle table

% ADD EVENT TO DATA
data.events=[StartofRec;EndofRec];
clear StartofRec EndofRec

% EXTRACT LFP FROM ALL TETRODES FOR LATER USE
% locate .mat file from previous run
disp('Loading LFP')
% computer=cd;
% computer=strsplit(computer,filesep);
% computer(2:end)=[];
% file=dir([computer{1},filesep,'**', filesep, data.rat,'_',data.sessionID,'.mat']);
% if ~isempty(file)
%     load([file.folder,filesep,file.name,],'lfp')
%     lfp.ts=lfp.ts.*10^6;
%     data.lfp=lfp;
%     clear lfp
% else
    channels=table2cell(struct2table(dir([path,filesep,'*.ncs'])));
    eegfile=[strcat(channels(:,2),filesep,channels(:,1))];
    [data]=downsampleLFP(eegfile,data);
    clear tetrodes eegfile
% end
data.lfp.lfpsamplerate=1000;

% FIX TIMESTAMPS (set a zero point and convert microseconds to seconds)
data=FixTime(data);

% START OF EVENT LOOP 1:nSESSIONS
for event=1:size(data.events,2)
    % SET UP SOME BASIC INFO ABOUT THE SESSION
    % LINEAR TRACK OR NOT
    if event==1 && isequal(linear_track,'yes')
        linear_track='yes';
    else
        linear_track = 'no';
    end
    
    % RESTRICT DATA BY START AND END OF EVENT TO CALCULATE SESSION
    % DURATION, AND MOVEMENT STATS
    data_video=data.frames;
    data_video=data_video(data_video(:,1)>=data.events(1,event)...
        & data_video(:,1)<=data.events(2,event),:);
    
    % SMOOTH VELOCITY WITH A WIDTH OF .8 SECONDS (median is use to limit outliers)
    data_video(:,5)=smoothdata(data_video(:,5),'movmedian',sampleRate*.8);

    % NEW LENGTH OF DATA
    lengthofdata=length(data_video);
    
    % DURATION OF SESSION (MIN)
    sessionduration=(lengthofdata/sampleRate)/60;
    data.session_duration(event)=sessionduration;
    disp(['SESSION_',num2str(event),' WAS ',num2str(sessionduration),' MIN'])
    if sessionduration<3
        disp('SESSION TOO SHORT...SKIPPING');
        continue;
    end
    clear sessionduration
    
    % Arena diameter
    if isequal(linear_track,'no') && event==4
        track_length = 100;
    elseif isequal(linear_track,'no') && sum(contains(path,'HPCatn'))==0
        track_length = 76.5;
    elseif isequal(linear_track,'no') && sum(contains(path,'HPCatn'))>0 && event>1
        track_length=100;
    end
    data.maze_size_cm(event)=track_length;
    
    % GET VELOCITY
    vel_cmPerSec=data_video(:,5);
        
    % EXTRACT BASIC MOVEMENT DATA
    data.BasicLoco.AverageAnglePerSec(event)=rad2deg(circ_mean((abs(diff(deg2rad(data_video(:,4)))))*sampleRate));
    data.BasicLoco.OverallDistTraveled(event)=nansum((vel_cmPerSec/sampleRate)*pixelDist);
    data.BasicLoco.MeanVelocity(event)=nanmean(vel_cmPerSec);

    % calculate running speed over 200ms bins
    n=round(.200*sampleRate);% average every n values
    data.binned_vel{event}=arrayfun(@(i) mean(vel_cmPerSec(i:i+n-1)),1:n:length(vel_cmPerSec)-n+1)';

    clear vel_cmPerSec vel_abs ExtractedAngle BasicLoco data_video n
    
    % START OF MAIN SPIKE TO PATH AND DATA ANALYSIS LOOP
    for i=1:length(data.Spikes)
        disp([extractBefore(data.spikesID.TetrodeNum{i},'.'),...
            ' Cell ',num2str(data.spikesID.CellNum(i))])
        
        [data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,event,i);
        
        % SEE IF SPIKES REGULARLY OCCUR OVER TIME
        temporal_stability=temporalstability(data_video_spk(data_video_spk(:,6)==1,1),data_video_spk(:,1));
        
        % calculate IFR
        ifr=IFR(data_video_spk,sampleRate);
        % correlation ifr to average running speed to see if firing rate
        % increases as a function of running speed
        b=data.binned_vel{event}; % the averaged vector
        % make sure two vectors are same length
        if length(ifr)~=length(b)
            [~,I]=max([length(ifr);length(b)]);
            if I==1
                ifr(length(b)+1:end)=[];
            elseif I==2
                b(length(ifr)+1:end)=[];
            end
        end
        SpeedScore=corr(ifr,b,'rows','complete');
        data.ifr{i,event}=ifr;
        clear ifr b I
        
        % SPLIT UP DATA BY RIGHT AND LEFT RUNS FOR LINEAR TRACK DATA
        if isequal(linear_track,'yes')
            runiteration=2;
            [right,left,DirectionalityIndex,Displacement,nlaps,rateoverlap,...
                fieldoverlap,spatialcorrelation,startgoodrun,stopgoodrun,laps]=...
                RightVsLeft(data_video_nospk,data_video_spk,track_length,sampleRate,data);
            if ~isfield(data.linear_track,'lapinfo')
                data.linear_track.lapinfo.startgoodrun=startgoodrun;
                data.linear_track.lapinfo.stopgoodrun=stopgoodrun;
                data.linear_track.lapinfo.laps=laps;
            end
            data.linear_track.right{i}.dataspks=right.dataspks;
            data.linear_track.right{i}.maps=right.maps;
            data.linear_track.right{i}.laps=right.laps;
            data.linear_track.left{i}.dataspks=left.dataspks;
            data.linear_track.left{i}.maps=left.maps;
            data.linear_track.left{i}.laps=left.laps;  
        else
            runiteration=1;
        end
        for iruns=1:runiteration
            if isequal(linear_track, 'yes') % UNPACK STRUCTURE
                if iruns==1
                    occ=right.occ; nBinsx=right.nBinsx; nBinsy=right.nBinsy;
                    Coherence=right.Coherence;
                    SmoothRateMap=right.SmoothRateMap;
                    
                    [fields]=getPlaceFields(SmoothRateMap,'minPeakRate',1,...
                        'minFieldWidth',3,'percentThreshold',.2,'maxFieldWidth',length(SmoothRateMap));

                    
                    data.linear_track.right{1, i}.fields=fields{1, 1};
                    data_video_smoothfilttemp=right.dataspks;
                    lap_perm_stability=right.lap_perm_stability;
                    stabilityoverlaps=right.stabilityoverlaps;
                    meanstability=right.meanstability;
                elseif iruns==2
                    occ=left.occ; nBinsx=left.nBinsx; nBinsy=left.nBinsy;
                    Coherence=left.Coherence;
                    SmoothRateMap=left.SmoothRateMap;
                    
                    [fields]=getPlaceFields(SmoothRateMap,'minPeakRate',1,...
                        'minFieldWidth',3,'percentThreshold',.2,'maxFieldWidth',length(SmoothRateMap));
                    
                    data.linear_track.left{1, i}.fields=fields{1, 1};
                    data_video_smoothfilttemp=left.dataspks;
                    lap_perm_stability=left.lap_perm_stability;
                    stabilityoverlaps=left.stabilityoverlaps;
                    meanstability=left.meanstability;
                end
                
                data_video_spk=data_video_smoothfilttemp;
                clear ia data_video_smoothfilttemp
            end
            spks_VEL=data_video_spk(data_video_spk(:,6)==1,:);            
            
            % THETA MODULATION
            [thetaindex,thetapeak,cor,~]=thetamodulation(spks_VEL(:,1));
            
            % BIN & SMOOTH OPEN FIELD DATA
            if isequal(linear_track,'no')
                [SmoothRateMap,nBinsx,nBinsy,occ,Coherence]=bindata(data_video_nospk,...
                    sampleRate,spks_VEL,linear_track,track_length);
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
                if event==3 && contains(data.mazetypes{event-1},'Cylinder') && contains(path,'PAE')
                    [Displacement,DisplacementCorr]=Displacement2(SmoothRateMap,...
                        data.ratemap{i,event});
                elseif event==3 && contains(data.mazetypes{event-1},'Cylinder')
                    [Displacement,DisplacementCorr]=Displacement2(SmoothRateMap,...
                        data.ratemap{i,event-1});
                end
                
                [borderScore,Field2Wall]=BorderScore(SmoothRateMap2);
                Field2Wall=Field2Wall*(track_length/length(SmoothRateMap));
                
                % find field
                [row,col]=find(field);k=boundary(row,col);
                bound=[col(k),row(k)];
                
                if length(bound)>2
                    E=Eccentricity(col(k),row(k));
                else
                    E=NaN;
                end
                
                % format data in a open field for phase precession                 
                [~, filename] = fileparts(ID{i});
                trodeID = str2double(extractAfter(filename,'TT'));
                
                % Implementation of the pass index technique for examination of open-field phase precession
                binside=mean([range(data_video_nospk(:,2))/length(SmoothRateMap),...
                    range(data_video_nospk(:,3))/length(SmoothRateMap)]);
                if sum(data_video_spk(:,6)==1)<5
                    occ4Ph=NaN;
                    fieldbound=[0 1];
                    spks_VEL4LFP=NaN;
                else
                    try
                    results=pass_index(data_video_nospk(:,1),data_video_nospk(:,2:3),...
                        data_video_spk(data_video_spk(:,6)==1,1),...
                        [data.lfp.ts(data.lfp.ts>=data.events(1,event) & data.lfp.ts<=data.events(2,event))]',...
                        [data.lfp.signal(trodeID,data.lfp.ts>=data.events(1,event) & data.lfp.ts<=data.events(2,event))]',...
                        'plots',0,'method','place','binside',binside);
                    occ4Ph=[results.ts,results.pass_index,zeros(length(results.ts),1)];
                    fieldbound=[0 1];
                    spks_VEL4LFP=data_video_spk(data_video_spk(:,6)==1,1);
                    catch
                        [spks_VEL4LFP,occ4Ph,fieldbound]=CylindarForPhasePrec(data_video_spk,SmoothRateMap,bound);
                    end
                end
                % calculate egocentric and border modulation based on (Peyrache, Schieferstein, Buzsaki, 2017)
                if contains(data.mazetypes{event},'Box')
                    out=egocentricmodulation(data_video_spk(:,2),data_video_spk(:,3),...
                        data_video_spk(:,4),data_video_spk,length(SmoothRateMap),sampleRate,3);
                    egomod=max(out);
                    
                    bordermod=bordermodulation(length(SmoothRateMap),data_video_spk(:,1),...
                        data_video_spk(:,2),data_video_spk(:,3),sampleRate,logical(data_video_spk(:,6)),3);
                    bordermod=max(bordermod);
                else
                    out=egocentricmodulation_circ(data_video_spk(:,2),data_video_spk(:,3),...
                    data_video_spk(:,4),data_video_spk,length(SmoothRateMap),data.samplerate,3);
                    egomod=max(out);
                    bordermod=bordermodulation_circ(SmoothRateMap,data_video_spk(:,1),...
                        data_video_spk(:,2),data_video_spk(:,3),data.samplerate,logical(data_video_spk(:,6)),3);
                    bordermod=max(bordermod);

                end
            else
                % LINEAR
                for fie=1:length(fields{1})
                    PR(fie)=fields{1}{fie}.peakFR;
                    start_(fie)=fields{1}{fie}.start;
                    stop_(fie)=fields{1}{fie}.stop;
                    temp=rescale([1:length(SmoothRateMap), fields{1, 1}{fie}.start,fields{1, 1}{fie}.stop],0,1);
                    fieldbound(fie,:)=temp(end-1:end);
                end

                [~,I]=max(PR);
                field=zeros(1,length(SmoothRateMap));
                field(fields{1, 1}{I}.start:fields{1, 1}{I}.stop)=1;
                FieldWidth=fields{1, 1}{I}.width*(track_length/length(field));
                
                clear temp PR

                % in field / out field FR
                Map4bordertemp=SmoothRateMap;Map4bordertemp(~field)=NaN;
                infieldFR=nanmean(reshape(Map4bordertemp,[],1));if isnan(infieldFR)==1;infieldFR=0;end
                Map4bordertemp=SmoothRateMap;Map4bordertemp(logical(field))=NaN;
                outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
                clear Map4bordertemp
                
                % spikes in field
                rescale_x=rescale(data_video_spk(:,2),0,1);
                nspikes_infield=sum(data_video_spk(rescale_x>fieldbound(1) & rescale_x<fieldbound(2),6));
            end
            
            % LFP ANALYSIS
            % Create name of CSC file from cell and path name
            [~, filename] = fileparts(ID{i});
            trodeID = str2double(extractAfter(filename,'TT'));

            if size(data.events,2)>1
                theta_phase=data.lfp.theta_phase(trodeID,...
                    data.lfp.ts>=data.events(1,event) &...
                    data.lfp.ts<=data.events(2,event));
                ts=data.lfp.ts(1,...
                    data.lfp.ts>=data.events(1,event) &...
                    data.lfp.ts<=data.events(2,event));
                theta=data.lfp.theta(trodeID,...
                    data.lfp.ts>=data.events(1,event) &...
                    data.lfp.ts<=data.events(2,event));                
            else
              theta_phase=data.lfp.theta_phase(trodeID,:);
              ts=data.lfp.ts;
              theta=data.lfp.theta(trodeID,:);
            end
            
            
            
            %   * POS_TS: Vector of time stamps for the sample state
            %   * POS: MXN matrix of the sample state, where M is the number of samples
            %   and N is the dimensions of POS
            %   * SPK_TS: Spike times for the cell
            %   * LFP_TS: Time stamps for the local field potential (LFP) Sample
            %   * LFP_SIG: The LFP signal
%             results = pass_index( varargin )
%             figure
%             binside=mean([range(data_video_nospk(:,2))/length(SmoothRateMap),...
%                 range(data_video_nospk(:,3))/length(SmoothRateMap)]);
%             
%             results=pass_index(data_video_nospk(:,1),[data_video_nospk(:,2),data_video_nospk(:,2)],...
%                 data_video_spk(data_video_spk(:,6)==1,1),...
%                 [data.lfp.ts(data.lfp.ts>=data.events(1,event) & data.lfp.ts<=data.events(2,event))]',...
%                 [data.lfp.signal(trodeID,data.lfp.ts>=data.events(1,event) & data.lfp.ts<=data.events(2,event))]',...
%                 'plots',1,'method','place','binside',binside);
            
            
            % calc phase precession stats
            if isequal(linear_track,'yes')
                % for linear track lets look at precession for every existing field
                for f=1:size(fieldbound)
                    rescale_x=rescale(data_video_spk(:,2),0,1);
                    temp=data_video_spk(rescale_x>fieldbound(f,1) & rescale_x<fieldbound(f,2),:);
                    [ThPrecess]=PHprecession([ts',theta_phase'],temp(temp(:,6)==1,:),temp(temp(:,6)==0,1:3),[0 1]);
                    if iruns==1
                        data.linear_track.right{1, i}.fields{f}.ThPrecess=ThPrecess;
                    elseif iruns==2
                        data.linear_track.left{1, i}.fields{f}.ThPrecess=ThPrecess;
                    end
                end
                % to retain our current measures stucture, lets find (variable I)
                % the phase precession output for the field with the highest peak
                if iruns==1
                    ThPrecess=data.linear_track.right{1, i}.fields{I}.ThPrecess;
                elseif iruns==2
                    ThPrecess=data.linear_track.left{1, i}.fields{I}.ThPrecess;
                end
            else
                [ThPrecess]=PHprecession([ts',theta_phase'],spks_VEL4LFP,occ4Ph,fieldbound);
            end
            % theta power
            
            [MeanThetaFreq,MeanOverallPow]=meanfreq(theta,1000);
           
            
            clear tetrode eegfile theta_phase ts fieldbound
            
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
            
            % GRID CELL ANALYSIS / GRID SCORE
            %             if isequal(linear_track,'no');
            %                 autocorr=SpatialAutoCorr(SmoothRateMap,length(SmoothRateMap));
            %                 gridout=GridScore_Sinusoidal(autocorr,length(SmoothRateMap));
            %                 gridscore=gridout.maxSinuGrid;
            %             end
            
            % INTRA-TRIAL STABILITY
            IntraTrialR=IntraTrialStability(data_video_spk,linear_track,track_length);
            
            % SPIKE DIRECTION
            [r,~,Direct_infoContent,~,preferred_Direction,~]=...
                tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),spks_VEL(:,4),sampleRate);
            
            % BURST INDEX
            [spkBurstIx,burstLg,burstIdx]=BurstSpikes(spks_VEL(:,1));
            
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
                borderScore=NaN; 
                E=NaN;
                DisplacementCorr=NaN;
                egomod=NaN;
                bordermod=NaN;
            end
            if ~exist('Displacement','var');Displacement=NaN;DisplacementCorr=NaN;end
            
            measures=[InformationContent,Coherence,sparsity,PeakRate,...
                OverallFR,Field2Wall,FieldWidth,nSpikes,r,...
                preferred_Direction(1),Direct_infoContent,DirectionalityIndex,...
                ThPrecess.slope,ThPrecess.RSquared,ThPrecess.lapSlope,ThPrecess.lapR2,...
                ThPrecess.lapCorrelation,ThPrecess.lapPhaseOffset,ThPrecess.circLinCorr,...
                ThPrecess.pval,ThPrecess.slopeCpU,ThPrecess.phaseOffset,...
                MeanThetaFreq,MeanOverallPow,...
                Displacement,infieldFR,outfieldFR,grade(i),borderScore,E,NumbActiveBins,...
                DisplacementCorr,SpeedScore,IntraTrialR,nspikes_infield,...
                ThPrecess.phaselock.Rlength,ThPrecess.phaselock.Pval,thetaindex,thetapeak,...
                clusterquality(i,:),prop(i,:),temporal_stability,nlaps,rateoverlap,fieldoverlap,...
                lap_perm_stability,stabilityoverlaps,meanstability,spatialcorrelation,egomod,bordermod,...
                burstIdx];
            
            if event==1 && iruns==1
                data.measures(i,:,event)=measures;
                data.ratemap{i,event}=SmoothRateMap;
                data.thetaautocorr{i,event}=cor;
                data.ThPrecess{i,event}=ThPrecess.scatteredPH;
            elseif event==1 && iruns>1
                data.measures(i,:,event+1)=measures;
                data.ratemap{i,event+1}=SmoothRateMap;
                data.thetaautocorr{i,event+1}=cor;
                data.ThPrecess{i,event+1}=ThPrecess.scatteredPH;
            elseif event>1 && ~contains(path,'ClarkP30_Recordings')
                data.measures(i,:,event+1)=measures;
                data.ratemap{i,event+1}=SmoothRateMap;
                data.thetaautocorr{i,event+1}=cor;
                data.ThPrecess{i,event+1}=ThPrecess.scatteredPH;
            elseif event>1 && contains(data.mazetypes{event-1},'Cylinder') 
                data.measures(i,:,event)=measures;
                data.ratemap{i,event}=SmoothRateMap;
                data.thetaautocorr{i,event}=cor;
                data.ThPrecess{i,event}=ThPrecess.scatteredPH;
            end
        end
        clearvars -except path data_video SpikeFile mclustpath tfile track_length...
            linear_track SmoothRateMap_Right SmoothRateMap_Left figures S...
            sampleRate i Overall_DirectionalityIndex event...
            timestamps vel_cmPerSec pixelDist grade data ratID video clusterquality prop ID...
            EEG_DownSampledTimestamps EEGthetaData EEG_DownSampledData root
    end
end
clearvars -except data figures 

% save mat file to temp file
if contains(data.session_path,'HPCatn')
    save(['D:\Projects\HPCatn\ProcessedData\',data.rat,'_',data.sessionID],'-struct','data','-v7.3')
elseif contains(data.session_path,'PAE')
    save(['D:\Projects\PAE_PlaceCell\ProcessedData\',data.rat,'_',data.sessionID],'-struct','data','-v7.3')
elseif contains(data.session_path,'ClarkP30')
    save(['F:\ClarkP30_Recordings\ProcessedData\',data.rat,'_',data.sessionID],'-struct','data','-v7.3')
end
% -------------------------------CREATING PLOTS----------------------------
if figures==1
    close all
    postprocessFigures(data);
end
disp 'DONE'

