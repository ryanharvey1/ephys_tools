% High level function that processes single unit and lfp data
%
% output variables listed in variable 'varnames' below
%
% Ryan E. Harvey
%
function data=postprocess(varargin)
p = inputParser;
p.addParameter('path',pwd);
p.addParameter('figures',1);
p.addParameter('manual_xy_cleanup',0);
p.addParameter('overwrite_lfp',0);
p.parse(varargin{:});

path = p.Results.path;
figures = p.Results.figures;
manual_xy_cleanup = p.Results.manual_xy_cleanup;
overwrite_lfp = p.Results.overwrite_lfp;


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

data.offset = ts(1);

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
data.basename = ratID{end};

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
data.mazetypes=get_maze_type(data);

% SET UP VARIABLE NAMES & DATA STRUCTURE
data.varnames={'InformationContent','Coherence','Sparsity','PeakRate',...
    'OverallFiringRate','Field2Wall','FieldWidth','nfields','nSpikes','mean_vector_length',...
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
    'egomod','bordermod','burstIdx',...
    'corrected_r','corrected_dic','corrected_info_content','corrected_sparsity','correction_fit'};

data.Spikes=S;
data.spikesID.paths=ID;
data.spikesID.TetrodeNum=tetrode;
data.spikesID.CellNum=cellnum;
data.avgwave=avgwave;
clear varnames S avgwave

data = get_session_idx(data);

% CHECK FOR MANUAL COORDINATE CORRECTIONS
if exist([path,filesep,'restrictxy.mat'],'file') || manual_xy_cleanup
    % run the following line if tracker errors from unplugs are present
    if manual_xy_cleanup
        manual_trackerjumps(ts,x,y,StartofRec,EndofRec,path);
    end
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

% make circular and linear track coords linear
if any(contains(data.mazetypes,'track','IgnoreCase',true))
    for nmaze=find(contains(data.mazetypes,'track','IgnoreCase',true))
        templinear=data_video(data_video(:,1)>=StartofRec(nmaze) & data_video(:,1)<=EndofRec(nmaze),:);
        data.linear_track{nmaze}.nonlinearFrames=templinear;
        if contains(data.mazetypes{nmaze},'circ track')
            [X,Y]=ParamToHorseshoeTrack(templinear(:,2),templinear(:,3));
            data_video(data_video(:,1)>=StartofRec(nmaze) & data_video(:,1)<=EndofRec(nmaze),2:3)=[X Y];
        elseif contains(data.mazetypes{nmaze},'Linear','IgnoreCase',true)
            [~,XY,~] = pca([templinear(:,2),templinear(:,3)]);
            data_video(data_video(:,1)>=StartofRec(nmaze) & data_video(:,1)<=EndofRec(nmaze),2:3)=XY;
        end
    end
end

% FIND MAX DIM OF MAZE TO CALCULATE VELOCITY
load('mazesize.mat','mazesize')
[~,b]=ismember(data.mazetypes,[mazesize{:,1}]');
data.maze_size_cm=[mazesize{b,2}];
clear b

% patch for linear track maze size bug
if datenum([data.sessionID(2:5),'-',data.sessionID(6:7),'-',data.sessionID(8:9)],'yyyy-mm-dd') <...
        datenum('2016-07-14','yyyy-mm-dd') && any(contains(data.mazetypes,'track','IgnoreCase',true))
    data.maze_size_cm(1) = 90;
end
mazesize=max(data.maze_size_cm);

% GET VELOCITY
[vel_cmPerSec,~,pixelDist]=InstaVel([data_video(:,2),data_video(:,3)],mazesize,data.samplerate);

% vel_abs=smooth(vel_abs,24);
data_video=[data_video,[vel_cmPerSec(1);vel_cmPerSec]];

data.frames=data_video;
clear data_video X Y vel_cmPerSec templinear XY x y ts angle StartofRec EndofRec

data = get_lfp(data,'overwrite_lfp',overwrite_lfp);

% FIX TIMESTAMPS (set a zero point and convert microseconds to seconds)
data=FixTime(data);

sess_idx = 1;
% START OF EVENT LOOP 1:nSESSIONS
for event=1:size(data.events,2)
    % SET UP SOME BASIC INFO ABOUT THE SESSION
    % LINEAR TRACK OR NOT
    if contains(data.mazetypes{event},'track','IgnoreCase',true)
        track=1;
    else
        track=0;
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
%     if data.session_duration(event)<3
%         disp('SESSION TOO SHORT...SKIPPING');
%         continue;
%     end
    
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
        [ifr,~]=instantfr(unique(data.Spikes{i}(data.Spikes{i}>data.events(1,event)...
            & data.Spikes{i}<data.events(2,event))),data_video_spk(1,1):.2:data_video_spk(end,1));
        
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
        if track==1
            
            runiteration=2;
            
            [right,left,DirectionalityIndex,Displacement,nlaps,rateoverlap,...
                fieldoverlap,spatialcorrelation,startgoodrun,stopgoodrun,laps]=...
                RightVsLeft(data_video_nospk,data_video_spk,data.maze_size_cm(event),data.samplerate,data,event);
            
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
            if track==1 % UNPACK STRUCTURE
                direction={'right','left'};
                
                occ=splitruns.(direction{iruns}).occ;
                nBinsx=splitruns.(direction{iruns}).nBinsx;
                nBinsy=splitruns.(direction{iruns}).nBinsy;
                ratemap=splitruns.(direction{iruns}).SmoothRateMap;
                Coherence = compute_spatial_coherence(ratemap);
                
                if sum(ratemap)==0
                    laps_ = ratemap';
                else
                    laps_=reshape([splitruns.(direction{iruns}).maps{:}],[],length(splitruns.(direction{iruns}).maps));
                end
                if isempty(laps_)
                    laps_ = ratemap';
                end
                [fields]=place_cell_analysis.getPlaceFields('ratemap',laps_,'minPeakRate',1,...
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
            if track==0
                [ratemap,nBinsx,nBinsy,occ,~]=bindata(data_video_nospk,...
                    data.samplerate,spks_VEL,track,data.maze_size_cm(event));
                SmoothRateMap2=ratemap;
                
                [Coherence] = compute_spatial_coherence(ratemap);
                
                fields = place_cell_analysis.getPlaceFields_2d('ratemap',ratemap,...
                    'maxFieldWidth',length(ratemap),...
                    'maze_size_cm',data.maze_size_cm(event),...
                    'debugging_fig',0);
                field = fields.masked_field{1};
                data.openfield{event}.fields{i} = fields;
                
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
                Field2Wall=Field2Wall*(data.maze_size_cm(event)/length(ratemap));
                
                % find field
                [row,col]=find(field);k=boundary(row,col);
                bound=[col(k),row(k)];
                
                if length(bound)>2
                    E=Eccentricity(col(k),row(k));
                else
                    E=NaN;
                end
                
                % format data in a open field for phase precession
                trodeID = get_channel_from_tetrode(data,data.spikesID.paths{i});
                
                % Implementation of the pass index technique for examination of open-field phase precession
                binside=mean([range(data_video_nospk(:,2))/length(ratemap),...
                    range(data_video_nospk(:,3))/length(ratemap)]);
                if sum(data_video_spk(:,6)==1)<10 || nansum(ratemap(:))==0
                    occ4Ph=NaN;
                    fieldbound=[0 1];
                    spks_VEL4LFP=NaN;
                else
                    results=pass_index(data_video_nospk(:,1),data_video_nospk(:,2:3),...
                        data_video_spk(data_video_spk(:,6)==1,1),...
                        [data.lfp.ts(data.lfp.ts>=data.events(1,event) & data.lfp.ts<=data.events(2,event))]',...
                        [data.lfp.signal(trodeID,data.lfp.ts>=data.events(1,event) & data.lfp.ts<=data.events(2,event))]',...
                        'plots',0,'method','place','binside',round(binside),'sample_along','arc_length');
                    
                    occ4Ph=[results.ts,results.pass_index,zeros(length(results.ts),1)];
                    fieldbound=[0 1];
                    spks_VEL4LFP=data_video_spk(data_video_spk(:,6)==1,1);
                end
                egomod=NaN;
                bordermod=NaN;
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
                FieldWidth=fields{1, 1}{I}.width*(data.maze_size_cm(event)/length(field));
                
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
            trodeID = get_channel_from_tetrode(data,data.spikesID.paths{i});

            if size(data.events,2)>1
                theta_phase=data.lfp.theta_phase(trodeID,...
                    data.lfp.ts>=data.events(1,event) &...
                    data.lfp.ts<=data.events(2,event));
                ts=data.lfp.ts(1,...
                    data.lfp.ts>=data.events(1,event) &...
                    data.lfp.ts<=data.events(2,event));
            else
                theta_phase=data.lfp.theta_phase(trodeID,:);
                ts=data.lfp.ts;
            end
            
            % calc phase precession stats
            if track==1
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
            [MeanThetaFreq,MeanOverallPow]=...
                meanfreq(data.lfp.signal(trodeID,data.lfp.ts>=data.events(1,event) &...
                data.lfp.ts<=data.events(2,event)),data.lfp.lfpsamplerate,[4,12]);
            
            clear tetrode eegfile theta_phase ts fieldbound
            
            InformationContent=place_cell_analysis.SpatialInformation('ratemap',...
                ratemap,'occupancy',occ,'n_spikes',sum(data_video_spk(:,6)));
            
            
            % Calculate the distance from place field to wall - peaks method
            if track==1
                [~, x_max]=find(ratemap==max(ratemap(:)));
                Field2Wall=min([nBinsx-x_max x_max-0])*(data.maze_size_cm(event)/length(ratemap));
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
            IntraTrialR=place_cell_analysis.IntraTrialStability(data_video_spk,track,data.maze_size_cm(event));
            
            % SPIKE DIRECTION
            [r,~,Direct_infoContent,~,preferred_Direction,hdTuning]=...
                tuningcurve(data_video_spk(data_video_spk(:,6)==0,4),spks_VEL(:,4),data.samplerate);
            
            % BURST INDEX
            [spkBurstIx,burstLg,burstIdx]=BurstSpikes(spks_VEL(:,1));
            
            
            % Correct for movement artifacts that may modulate position or
            % direction related firing using factorial Maximum Likelihood Model
            if track==0
                
                [corrected_p,corrected_d,correction_fit]=FMLM_wrapper(data,i,event);
                
                if ~isnan(correction_fit)
                    
                    % corrected direction metrics
                    corrected_r = circ_r(deg2rad(3:6:357)',corrected_d,deg2rad(6));
                    corrected_dic = HD_cell_analysis.computeDIC(histcounts(data_video_spk(data_video_spk(:,6)==0,4),0:6:360),...
                        corrected_d',sum(data_video_spk(:,6))/sum(data_video_spk(:,6)==0)*data.samplerate);
                    
                    % corrected position metrics
                    corrected_info_content = place_cell_analysis.SpatialInformation('ratemap',...
                        corrected_p,'occupancy',occ,'n_spikes',sum(data_video_spk(:,6)));
                    corrected_sparsity = place_cell_analysis.Sparsity('ratemap',corrected_p,'occupancy',occ);
                else
                    corrected_r=NaN;
                    corrected_dic=NaN;
                    corrected_info_content=NaN;
                    corrected_sparsity=NaN;
                    correction_fit=NaN;
                end
            end
            % pack results for NaN measures given maze type
            if track==0
                DirectionalityIndex=NaN;
                nlaps=NaN;
                rateoverlap=NaN;
                fieldoverlap=NaN;
                lap_perm_stability=NaN;
                stabilityoverlaps=NaN;
                meanstability=NaN;
                spatialcorrelation=NaN;
                
                FieldWidth = fields.fieldwidth{1};
                nfields = fields.nfields;
            else
                borderScore=NaN;
                E=NaN;
                DisplacementCorr=NaN;
                egomod=NaN;
                bordermod=NaN;
                corrected_r=NaN;
                corrected_dic=NaN;
                corrected_info_content=NaN;
                corrected_sparsity=NaN;
                correction_fit=NaN;
                corrected_d=NaN;
                nfields = length(fields{1,1});
            end
            if ~exist('Displacement','var');Displacement=NaN;DisplacementCorr=NaN;end
            
            % compile measures into vector (would be better as a table in the future)
            measures=[InformationContent,Coherence,sparsity,max(ratemap(:)),...
                nanmean(ratemap(:)),Field2Wall,FieldWidth,nfields,sum(data_video_spk(:,6)),r,...
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
                burstIdx,corrected_r,corrected_dic,corrected_info_content,corrected_sparsity,correction_fit];
            

            data.measures(i,:,data.session_idx(sess_idx))=measures;
            data.ratemap{i,data.session_idx(sess_idx)}=ratemap;
            data.thetaautocorr{i,data.session_idx(sess_idx)}=cor;
            data.ThPrecess{i,data.session_idx(sess_idx)}=ThPrecess.scatteredPH;
            data.hdTuning{i,data.session_idx(sess_idx)}=hdTuning;
            data.hdTuning_corrected{i,data.session_idx(sess_idx)}=corrected_d;
            
            sess_idx = sess_idx + 1;
        end
        clearvars -except track figures i event data clusterquality prop pixelDist path sess_idx
    end
end
clearvars -except data figures

data = orderfields(data);

% save mat file to processed data folder
processedpath=strsplit(data.session_path,filesep);
processedpath(end-2:end)=[];
save(fullfile(strjoin(processedpath,filesep),'ProcessedData',[data.rat,'_',data.sessionID]),'-struct','data','-v7.3')

% -------------------------------CREATING PLOTS----------------------------
if figures==1
    close all
    postprocessFigures.main(data,'colorcode','HD');
end
disp 'DONE'