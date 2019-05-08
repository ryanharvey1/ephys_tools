% High level function that processes single unit and lfp data
%
% output variables listed in variable 'varnames' below
%
% Ryan E. Harvey
%
function data=postprocess(path,mazesize,linear_track,figures)

addpath(fullfile('Analysis'),...
    genpath(fullfile('io')),...
    fullfile('LFP'),...
    genpath(fullfile('preprocessing')),...
    fullfile('Utils'),...
    fullfile('Visualize'),...
    fullfile('external_packages','CircStat2012a'),...
    genpath(fullfile('external_packages','FMAToolbox')),...
    fullfile('external_packages','Pass_Index'),...
    fullfile('external_packages','plotSpikeRaster'),...
    fullfile('external_packages','panel'),...
    fullfile('external_packages','Colormaps'),...
    fullfile('external_packages','Inpaint_nans'))


% load spike data
[S,avgwave,ID,cellnum,tetrode,clusterquality]=load_spikes(path);

disp(path)
if ~exist('S','var')
    disp('No Clusters... Skipping')
    return
end
disp(['    ',num2str(length(S)),' Cells'])

% CALCULATE WAVEFORM PROPERTIES
prop=waveformprop(avgwave);

% READ VIDEO DATA
[ts,x,y,angles] = Nlx2MatVT([path,filesep,'VT1.nvt'],[1,1,1,1,0,0],0,1);

% calculate sample rate from ts
data.samplerate=get_framerate(ts);

% DURATION OF SESSION (MIN)
if ((length(ts)/data.samplerate)/60)<4
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
% ADD EVENT TO DATA
data.events=[StartofRec;EndofRec];

% MAZE TYPES
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pull maze type from meta data if it exists
try
    data.mazetypes=get_maze_type(data);
catch
    for mt=1:length(StartofRec)
        if mt<5 && contains(path,'ClarkP30')
            mazetype{mt}='Cylinder';
        elseif mt>=5 && contains(path,'ClarkP30')
            mazetype{mt}=input(['Maze Type for Session: ',num2str(mt),'  '],'s');
        end
        if contains(linear_track,'yes') && ~contains(path,'HPCatn') && mt==1
            mazetype{mt}='LinearTrack';
        elseif contains(linear_track,'yes') && contains(path,'HPCatn') && mt==1
            mazetype{mt}='circ track';
        elseif contains(linear_track,'yes') && ~contains(path,'HPCatn') && mt>1
            mazetype{mt}='Cylinder';
        elseif contains(linear_track,'yes') && contains(path,'HPCatn') && mt>1
            mazetype{mt}='Box';
        end
    end
    data.mazetypes=mazetype; clear mt mazetypes
end

% SET UP VARIABLE NAMES & DATA STRUCTURE
data.varnames={'InformationContent','Coherence','Sparsity','PeakRate',...
    'OverallFiringRate','Field2Wall','FieldWidth','nSpikes','mean_vector_length',...
    'preferred_Direction','Direct_infoContent','DirectionalityIndex',...
    'PhPrecessSlope','PH_RSquared','lapPhPrecessSlope','lapPH_RSquared',...
    'PhlapCorrelation','lapPhaseOffset','PhcircLinCorr','Phpval','slopeCpU',...
    'phaseOffset','MeanThetaFreq','MeanThetaPow','Displacement',...
    'infieldFR','outfieldFR','borderScore','E','NumbActiveBins',...
    'DisplacementCorr','SpeedScore','IntraTrialStability','nspikes infield',...
    'ThetaPhaseLock R','ThetaPhaselock Pval','thetaindex','thetapeak',...
    'LRatio','ShortISI','IsolationDistance','NumberOfSpikes',...
    'WaveformPeak','WaveformValley','spikewidth','halfwidth','peak2valleyratio','peak2valleyslope','SpikeAmplitude',...
    'temporalstability','nlaps','rateoverlap','fieldoverlap',...
    'lap_perm_stability','stabilityoverlaps','meanstability','spatialcorrelation',...
    'egomod','bordermod','burstIdx'};

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
[xtemp,ytemp]=FixPos(x',y',ts',round(0.1667*data.samplerate));
tempangle=wrapTo360(fixNLXangle(angles',round(0.1667*data.samplerate)))'-90;
tempangle(tempangle<0)=tempangle(tempangle<0)+360;

data_video=[ts',xtemp,ytemp,tempangle]; clear xtemp ytemp ts

if isequal(linear_track,'yes')
    for nmaze=find(contains(data.mazetypes,'track','IgnoreCase',true))
        templinear=data_video(data_video(:,1)>=StartofRec(nmaze) & data_video(:,1)<=EndofRec(nmaze),:);
        data.linear_track{nmaze}.nonlinearFrames=templinear;
        if contains(data.mazetypes{nmaze},'circ track')
            [X,Y]=ParamToHorseshoeTrack(templinear(:,2),templinear(:,3));
        elseif contains(data.mazetypes{nmaze},'Linear','IgnoreCase',true)
            [X,Y]=ParameterizeToLinearTrack2(templinear(:,2),templinear(:,3));
        end
        data_video(data_video(:,1)>=StartofRec(nmaze) & data_video(:,1)<=EndofRec(nmaze),2:3)=[X Y];
    end
end

% FIND MAX DIM OF MAZE TO CALCULATE VELOCITY

if sum(contains(path,'HPCatn'))>0 && any(contains(data.mazetypes,'track'))
    mazesize=360;
elseif sum(contains(path,'HPCatn'))==0 && size(StartofRec,2)==1 && mazesize~=90
    mazesize=120;
elseif any(contains(data.mazetypes,'box'))
    mazesize=100;
elseif contains(linear_track,'no')
    mazesize=76.5;
end

% GET VELOCITY
[vel_cmPerSec,~,pixelDist]=InstaVel([data_video(:,2),data_video(:,3)],linear_track,mazesize,data.samplerate);

% vel_abs=smooth(vel_abs,24);
data_video=[data_video,[vel_cmPerSec(1);vel_cmPerSec]];

data.frames=data_video;
clear data_video X Y vel_cmPerSec templinear XY x y ts angle table


clear StartofRec EndofRec

% EXTRACT LFP FROM ALL TETRODES FOR LATER USE
disp('Have patience...Loading LFP')
channels=table2cell(struct2table(dir([path,filesep,'*.ncs'])));
eegfile=[strcat(channels(:,2),filesep,channels(:,1))];
[data]=downsampleLFP(eegfile,data);
clear tetrodes eegfile
data.lfp.lfpsamplerate=1000;

% FIX TIMESTAMPS (set a zero point and convert microseconds to seconds)

data=FixTime(data);

% START OF EVENT LOOP 1:nSESSIONS
for event=1:size(data.events,2)
    % SET UP SOME BASIC INFO ABOUT THE SESSION
    % LINEAR TRACK OR NOT
    if contains(data.mazetypes{event},'track','IgnoreCase',true)
        linear_track='yes';
    else
        linear_track='no';
    end
    
    % RESTRICT DATA BY START AND END OF EVENT TO CALCULATE SESSION
    % DURATION, AND MOVEMENT STATS
    data_video=data.frames(data.frames(:,1)>=data.events(1,event)...
        & data.frames(:,1)<=data.events(2,event),:);
    
    % SMOOTH VELOCITY WITH A WIDTH OF .8 SECONDS (median is use to limit outliers)
    data_video(:,5)=smoothdata(data_video(:,5),'movmedian',data.samplerate*.8);
    
    % DURATION OF SESSION (MIN)
    data.session_duration(event)=(length(data_video)/data.samplerate)/60;
    disp(['SESSION_',num2str(event),' WAS ',num2str(data.session_duration(event)),' MIN'])
    if data.session_duration(event)<3
        disp('SESSION TOO SHORT...SKIPPING');
        continue;
    end
    
    % Arena diameter
    if contains(data.mazetypes{event},'box','IgnoreCase',true)
        mazesize = 100;
    elseif contains(data.mazetypes{event},'Cylinder','IgnoreCase',true)
        mazesize = 76.5;
    end
    data.maze_size_cm(event)=mazesize;
    
    % GET VELOCITY
    vel_cmPerSec=data_video(:,5);
    
    % EXTRACT BASIC MOVEMENT DATA
    data.BasicLoco.AverageAnglePerSec(event)=rad2deg(circ_mean((abs(diff(deg2rad(data_video(:,4)))))*data.samplerate));
    data.BasicLoco.OverallDistTraveled(event)=nansum((vel_cmPerSec/data.samplerate)*pixelDist);
    data.BasicLoco.MeanVelocity(event)=nanmean(vel_cmPerSec);
    
    % calculate running speed over 200ms bins
    data.binned_vel{event}=binned_vel(data_video(:,1),vel_cmPerSec,.2);
    
    clear vel_cmPerSec vel_abs ExtractedAngle
    
    % START OF MAIN SPIKE TO PATH AND DATA ANALYSIS LOOP
    for i=1:length(data.Spikes)
        disp([extractBefore(data.spikesID.TetrodeNum{i},'.'),...
            ' Cell ',num2str(data.spikesID.CellNum(i))])
        
        % combine and limit video and spike data 
        [data_video_spk,data_video_nospk]=createframes_w_spikebinary(data,event,i);
        
        % SEE IF SPIKES REGULARLY OCCUR OVER TIME
        temporal_stability=temporalstability(data_video_spk(data_video_spk(:,6)==1,1),data_video_spk(:,1));
        
        % calculate IFR
        [ifr,~]=instantfr(data.Spikes{i}(data.Spikes{i}>data.events(1,event)...
            & data.Spikes{i}<data.events(2,event)),data_video_spk(1,1):.2:data_video_spk(end,1));
        
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
            data.binned_vel{event}=b;
        end
        SpeedScore=corr(b,ifr,'rows','complete');

        data.ifr{i,event}=ifr;
        clear ifr b I
        
        % SPLIT UP DATA BY RIGHT AND LEFT RUNS FOR LINEAR TRACK DATA
        if isequal(linear_track,'yes')
            
            runiteration=2;
            
            [right,left,DirectionalityIndex,Displacement,nlaps,rateoverlap,...
                fieldoverlap,spatialcorrelation,startgoodrun,stopgoodrun,laps]=...
                RightVsLeft(data_video_nospk,data_video_spk,mazesize,data.samplerate,data,event);
            
            if ~isfield(data.linear_track{event},'lapinfo')
                data.linear_track{event}.lapinfo.startgoodrun=startgoodrun;
                data.linear_track{event}.lapinfo.stopgoodrun=stopgoodrun;
                data.linear_track{event}.lapinfo.laps=laps;
            end
            data.linear_track{event}.right{i}.dataspks=right.dataspks;
            data.linear_track{event}.right{i}.maps=right.maps;
            data.linear_track{event}.right{i}.laps=right.laps;
            data.linear_track{event}.left{i}.dataspks=left.dataspks;
            data.linear_track{event}.left{i}.maps=left.maps;
            data.linear_track{event}.left{i}.laps=left.laps;
            
            splitruns.right=right;
            splitruns.left=left;

        else
            runiteration=1;
        end
        for iruns=1:runiteration
            if isequal(linear_track, 'yes') % UNPACK STRUCTURE
                direction={'right','left'};
                
                occ=splitruns.(direction{iruns}).occ; 
                nBinsx=splitruns.(direction{iruns}).nBinsx; 
                nBinsy=splitruns.(direction{iruns}).nBinsy;
                Coherence=splitruns.(direction{iruns}).Coherence;
                ratemap=splitruns.(direction{iruns}).SmoothRateMap;
                
                [fields]=place_cell_analysis.getPlaceFields(ratemap,'minPeakRate',1,...
                    'minFieldWidth',3,'percentThreshold',.2,'maxFieldWidth',length(ratemap));
                
                data.linear_track{event}.(direction{iruns}){1, i}.fields=fields{1, 1};
                data_video_smoothfilttemp=splitruns.(direction{iruns}).dataspks;
                lap_perm_stability=splitruns.(direction{iruns}).lap_perm_stability;
                stabilityoverlaps=splitruns.(direction{iruns}).stabilityoverlaps;
                meanstability=splitruns.(direction{iruns}).meanstability;
                
                data_video_spk=data_video_smoothfilttemp;
                clear ia data_video_smoothfilttemp
            end
            spks_VEL=data_video_spk(data_video_spk(:,6)==1,:);
            
            % THETA MODULATION
            [thetaindex,thetapeak,cor,~]=thetamodulation(spks_VEL(:,1));
            
            % BIN & SMOOTH OPEN FIELD DATA
            if isequal(linear_track,'no')
                [ratemap,nBinsx,nBinsy,occ,Coherence]=bindata(data_video_nospk,...
                    data.samplerate,spks_VEL,linear_track,mazesize);
                SmoothRateMap2=ratemap;
                [field,FieldWidth]=FindFF2D(SmoothRateMap2);
                FieldWidth=FieldWidth*(mazesize/length(field));
                if isnan(field)==0
                    % in field / out field FR
                    nanPlacement=isnan(SmoothRateMap2);
                    SmoothRateMap2(~field)=0;
                    SmoothRateMap2(nanPlacement)=NaN;
                    Map4bordertemp=SmoothRateMap2;
                    Map4bordertemp(~field)=NaN;
                    infieldFR=nanmean(reshape(Map4bordertemp,[],1));
                    Map4bordertemp=ratemap;
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
                    [Displacement,DisplacementCorr]=Displacement2(ratemap,...
                        data.ratemap{i,event});
                elseif event==3 && contains(data.mazetypes{event-1},'Cylinder')
                    [Displacement,DisplacementCorr]=Displacement2(ratemap,...
                        data.ratemap{i,event-1});
                end
                
                [borderScore,Field2Wall]=BorderScore(SmoothRateMap2);
                Field2Wall=Field2Wall*(mazesize/length(ratemap));
                
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
                binside=mean([range(data_video_nospk(:,2))/length(ratemap),...
                    range(data_video_nospk(:,3))/length(ratemap)]);
                if sum(data_video_spk(:,6)==1)<5
                    occ4Ph=NaN;
                    fieldbound=[0 1];
                    spks_VEL4LFP=NaN;
                else
                    results=pass_index(data_video_nospk(:,1),data_video_nospk(:,2:3),...
                        data_video_spk(data_video_spk(:,6)==1,1),...
                        [data.lfp.ts(data.lfp.ts>=data.events(1,event) & data.lfp.ts<=data.events(2,event))]',...
                        [data.lfp.signal(trodeID,data.lfp.ts>=data.events(1,event) & data.lfp.ts<=data.events(2,event))]',...
                        'plots',1,'method','place','binside',round(binside));
                    occ4Ph=[results.ts,results.pass_index,zeros(length(results.ts),1)];
                    fieldbound=[0 1];
                    spks_VEL4LFP=data_video_spk(data_video_spk(:,6)==1,1);
                end
                
                % calculate egocentric and border modulation based on (Peyrache, Schieferstein, Buzsaki, 2017)
                if contains(data.mazetypes{event},'Box')
                    out=egocentricmodulation(data_video_spk(:,2),data_video_spk(:,3),...
                        data_video_spk(:,4),data_video_spk,length(ratemap),data.samplerate,3);
                    egomod=max(out);
                    
                    bordermod=bordermodulation(length(ratemap),data_video_spk(:,1),...
                        data_video_spk(:,2),data_video_spk(:,3),data.samplerate,logical(data_video_spk(:,6)),3);
                    bordermod=max(bordermod);
                else
                    out=egocentricmodulation_circ(data_video_spk(:,2),data_video_spk(:,3),...
                        data_video_spk(:,4),data_video_spk,length(ratemap),data.samplerate,3);
                    egomod=max(out);
                    bordermod=bordermodulation_circ(ratemap,data_video_spk(:,1),...
                        data_video_spk(:,2),data_video_spk(:,3),data.samplerate,logical(data_video_spk(:,6)),3);
                    bordermod=max(bordermod);
                    
                end
            else
                % LINEAR
                for fie=1:length(fields{1})
                    PR(fie)=fields{1}{fie}.peakFR;
                    start_(fie)=fields{1}{fie}.start;
                    stop_(fie)=fields{1}{fie}.stop;
                    temp=rescale([1:length(ratemap), fields{1, 1}{fie}.start,fields{1, 1}{fie}.stop],0,1);
                    fieldbound(fie,:)=temp(end-1:end);
                end
                
                [~,I]=max(PR);
                field=zeros(1,length(ratemap));
                field(fields{1, 1}{I}.start:fields{1, 1}{I}.stop)=1;
                FieldWidth=fields{1, 1}{I}.width*(mazesize/length(field));
                
                clear temp PR
                
                % in field / out field FR
                Map4bordertemp=ratemap;Map4bordertemp(~field)=NaN;
                infieldFR=nanmean(reshape(Map4bordertemp,[],1));if isnan(infieldFR)==1;infieldFR=0;end
                Map4bordertemp=ratemap;Map4bordertemp(logical(field))=NaN;
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
            
            % calc phase precession stats
            if isequal(linear_track,'yes')
                % for linear track lets look at precession for every existing field
                for f=1:size(fieldbound)
                    rescale_x=rescale(data_video_spk(:,2),0,1);
                    
                    temp=data_video_spk(rescale_x>fieldbound(f,1) & rescale_x<fieldbound(f,2),:);
                    
                    [ThPrecess]=place_cell_analysis.PHprecession([ts',theta_phase'],temp(temp(:,6)==1,:),temp(temp(:,6)==0,1:3),[0 1]);
                    
                    data.linear_track{event}.(direction{iruns}){1, i}.fields{f}.ThPrecess=ThPrecess;
                end
                % to retain our current measures stucture, lets find (variable I)
                % the phase precession output for the field with the highest peak
                ThPrecess=data.linear_track{event}.(direction{iruns}){1, i}.fields{I}.ThPrecess;

            else
                [ThPrecess]=place_cell_analysis.PHprecession([ts',theta_phase'],spks_VEL4LFP,occ4Ph,fieldbound);
            end
            
            % theta power and freq
            [MeanThetaFreq,MeanOverallPow]=meanfreq(theta,1000);
            
            clear tetrode eegfile theta_phase ts fieldbound
            
            InformationContent=place_cell_analysis.SpatialInformation('ratemap',...
                ratemap,'occupancy',occ,'n_spikes',sum(data_video_spk(:,6)));
            
            
            % Calculate the distance from place field to wall - peaks method
            if isequal(linear_track,'yes')
                [~, x_max]=find(ratemap==max(ratemap(:)));
                Field2Wall=min([nBinsx-x_max x_max-0])*(mazesize/length(ratemap));
            end
            
            % Get Sparsity
            sparsity = place_cell_analysis.Sparsity('ratemap',ratemap,'occupancy',occ);
            
            
            % calculate percent of active bins (estimate of field size - Royer et al., 2010)
            NumbActiveBins=numel(find(ratemap > max(ratemap(:))*0.20));
            clear R rY r_max RC log_relR nBins relR occRSHP occSUM pX nCells x_max
            
            % GRID CELL ANALYSIS / GRID SCORE
            %             if isequal(linear_track,'no');
            %                 autocorr=SpatialAutoCorr(SmoothRateMap,length(SmoothRateMap));
            %                 gridout=GridScore_Sinusoidal(autocorr,length(SmoothRateMap));
            %                 gridscore=gridout.maxSinuGrid;
            %             end
            
            % INTRA-TRIAL STABILITY
            IntraTrialR=place_cell_analysis.IntraTrialStability(data_video_spk,linear_track,mazesize);
            
            % SPIKE DIRECTION
            [r,~,Direct_infoContent,~,preferred_Direction,~]=...
                tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),spks_VEL(:,4),data.samplerate);
            
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
            
            measures=[InformationContent,Coherence,sparsity,max(ratemap(:)),...
                nanmean(ratemap(:)),Field2Wall,FieldWidth,sum(data_video_spk(:,6)),r,...
                preferred_Direction(1),Direct_infoContent,DirectionalityIndex,...
                ThPrecess.slope,ThPrecess.RSquared,ThPrecess.lapSlope,ThPrecess.lapR2,...
                ThPrecess.lapCorrelation,ThPrecess.lapPhaseOffset,ThPrecess.circLinCorr,...
                ThPrecess.pval,ThPrecess.slopeCpU,ThPrecess.phaseOffset,...
                MeanThetaFreq,MeanOverallPow,...
                Displacement,infieldFR,outfieldFR,borderScore,E,NumbActiveBins,...
                DisplacementCorr,SpeedScore,IntraTrialR,nspikes_infield,...
                ThPrecess.phaselock.Rlength,ThPrecess.phaselock.Pval,thetaindex,thetapeak,...
                clusterquality(i,:),sum(data_video_spk(:,6)),prop(i,:),temporal_stability,nlaps,rateoverlap,fieldoverlap,...
                lap_perm_stability,stabilityoverlaps,meanstability,spatialcorrelation,egomod,bordermod,...
                burstIdx];
            
            if event==1 && iruns==1
                data.measures(i,:,event)=measures;
                data.ratemap{i,event}=ratemap;
                data.thetaautocorr{i,event}=cor;
                data.ThPrecess{i,event}=ThPrecess.scatteredPH;
            elseif event==1 && iruns>1
                data.measures(i,:,event+1)=measures;
                data.ratemap{i,event+1}=ratemap;
                data.thetaautocorr{i,event+1}=cor;
                data.ThPrecess{i,event+1}=ThPrecess.scatteredPH;
            elseif event>1 && contains(data.mazetypes{event},'track','IgnoreCase',true) && iruns==1
                data.measures(i,:,event+1)=measures;
                data.ratemap{i,event+1}=ratemap;
                data.thetaautocorr{i,event+1}=cor;
                data.ThPrecess{i,event+1}=ThPrecess.scatteredPH;
            elseif event>1 && contains(data.mazetypes{event},'track','IgnoreCase',true) && iruns==2
                data.measures(i,:,event+2)=measures;
                data.ratemap{i,event+2}=ratemap;
                data.thetaautocorr{i,event+2}=cor;
                data.ThPrecess{i,event+2}=ThPrecess.scatteredPH;
            elseif event>1 && contains(data.mazetypes{event},'Cylinder','IgnoreCase',true) && ~contains(data.session_path,'PAE')
                data.measures(i,:,event)=measures;
                data.ratemap{i,event}=ratemap;
                data.thetaautocorr{i,event}=cor;
                data.ThPrecess{i,event}=ThPrecess.scatteredPH;
            elseif event>1 && contains(data.mazetypes{event},'Cylinder','IgnoreCase',true) && contains(data.session_path,'PAE')
                data.measures(i,:,event+1)=measures;
                data.ratemap{i,event+1}=ratemap;
                data.thetaautocorr{i,event+1}=cor;
                data.ThPrecess{i,event+1}=ThPrecess.scatteredPH;
            elseif event>1 && contains(data.mazetypes{event},'Box','IgnoreCase',true) && ~contains(data.mazetypes{event-1},'track','IgnoreCase',true) 
                data.measures(i,:,event)=measures;
                data.ratemap{i,event}=ratemap;
                data.thetaautocorr{i,event}=cor;
                data.ThPrecess{i,event}=ThPrecess.scatteredPH;
            elseif event>1 && contains(data.mazetypes{event},'Box','IgnoreCase',true)
                data.measures(i,:,event+1)=measures;
                data.ratemap{i,event+1}=ratemap;
                data.thetaautocorr{i,event+1}=cor;
                data.ThPrecess{i,event+1}=ThPrecess.scatteredPH;
            end
        end
        clearvars -except path data_video SpikeFile mclustpath tfile mazesize...
            linear_track SmoothRateMap_Right SmoothRateMap_Left figures S...
            i Overall_DirectionalityIndex event...
            timestamps vel_cmPerSec pixelDist grade data ratID video clusterquality prop ID...
            EEG_DownSampledTimestamps EEGthetaData EEG_DownSampledData root
    end
end
clearvars -except data figures

% save mat file to processed data folder
processedpath=strsplit(data.session_path,filesep);
processedpath(end-2:end)=[];
save(fullfile(strjoin(processedpath,filesep),'ProcessedData',[data.rat,'_',data.sessionID]),'-struct','data','-v7.3')

% -------------------------------CREATING PLOTS----------------------------
if figures==1
    close all
    postprocessFigures.main(data);
end
disp 'DONE'

