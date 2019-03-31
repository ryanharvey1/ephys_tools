
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


function postprocess_snapsort_HD(path,track_length,linear_track,figures)
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

ns=1;grade=[];
for m=1:length(file)
    load(file{m});load(info{m})
    grade=[grade;final_grades'];
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

% READ VIDEO FILE
% Creates video tracking text file if one doesn't exist
if exist([path,filesep,'VT1.txt'],'file')==0; ReadVideoTrackingFile3(path); end

% locates newly created VT1.txt file
fileID=fopen([path,filesep,'VT1.txt'],'r');
dataArray=textscan(fileID,'%f%f%f%f%f%[^\n\r]','Delimiter','\t','TextType','string','EmptyValue',NaN,'ReturnOnError',false);
fclose(fileID);
video=[dataArray{1:end-1}];

% SPLIT UP SESSIONS BY EVENT FILES
[ timestamps,StartofRec,EndofRec ]=EventSplit(path);if isempty(timestamps);clear timestamps;end

% SET UP VARIABLE NAMES & DATA STRUCTURE
varnames={'InformationContent','Coherence','Sparsity','PeakRate',...
    'OverallFiringRate','Field2Wall','FieldWidth','nSpikes','mean_vector_length',...
    'preferred_Direction','Direct_infoContent','DirectionalityIndex',...
    'PhPrecessCorr','DepthofModulation','PhPrecessSlope','PhPrecessMeanFR','PH_RSquared',...
    'MeanThetaFreq','MeanOverallPow','ThetaRatio','PerSpkslessthan2ms','Displacement',...
    'infieldFR','outfieldFR','Cluster Grade','borderScore','E','NumbActiveBins',...
    'DisplacementCorr','SpeedScore','IntraTrialStability','nspikes infield',...
    'ThetaPhaseLock R','ThetaPhaselock Pval','thetaindex','thetapeak'};

ratID=strsplit(path,filesep);
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).varnames=varnames;
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).Spikes=cellfun(@(x) x*1000000,S,'un',0); %cell array of timestamps
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).avgwave=avgwave; %from snap

% ADD ENTIRE VIDEO TO DATA
[table,SumofNaN1,SumofNaN2,~]=NonDetects(video,length(video));
if SumofNaN1>SumofNaN2; VelIndex=[4 5]; else; VelIndex=[2 3]; end
[vel_cmPerSec,~,~] = InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,30);
table = table((vel_cmPerSec<=100),:);
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
    disp(['SESSION_',num2str(event),' WAS ',num2str(sessionduration),' MIN'])
    if sessionduration<2;disp('SESSION TOO SHORT...SKIPPING');continue;end
    
    % DEAL WITH NON DETECTS
    [table,SumofNaN1,SumofNaN2,~]=NonDetects(table,lengthofdata);
    
    % VELOCITY FILTER & SMOOTHING 
    % Arena diameter - INCLUDE LATER LB 11-Feb-2018 (see code graveyard -
    % may be helpful for radial arm maze)
    
    % CHOOSE LED WITH BEST TRACKING TO VELOCITY FILTER BY
    if SumofNaN1>SumofNaN2; VelIndex=[4 5]; else; VelIndex=[2 3]; end
    
    % CALCULATE VEL TO REMOVE JUMPS IN DATA
    [vel_cmPerSec,~,~] = InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,sampleRate);
    
    % find xy coords and spikes when rat is within velocity threshold ( <100)
    table = table((vel_cmPerSec<=100),:);
    
    filteredlengthofdata=length(table);
    disp(['filtering jumps in tracker ', num2str(lengthofdata-filteredlengthofdata),...
        ' data points, which is ', num2str(100-((filteredlengthofdata/lengthofdata)*100)), ' percent of your data.']);
    
    % SMOOTHING RAW XY DATA FROM EACH LED
    padDsData=[repmat(table(1,:),30,1); table; repmat(table(end,:),30,1)]; %pad ends with last data point
    padDsData=[padDsData(:,1),smooth(padDsData(:,2),10),smooth(padDsData(:,3),10),smooth(padDsData(:,4),10),smooth(padDsData(:,5),10)];
    table=(padDsData(30+1:end-30,:)); %remove pad
    
    % GET VELOCITY
    [vel_cmPerSec,vel_abs,pixelDist]=InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,sampleRate);
    
    % CONCAT TIMESTAMP/XY POSITION/ANGLE TOGETHER
    data_video=double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2)]);
    
    % STRAIGHTEN LINEAR TRACK
%     if isequal(linear_track,'yes')
%         [X,Y]=ParameterizeToLinearTrack2(data_video(:,2),data_video(:,3));data_video=[data_video(:,1) X Y];
%     end
    
    % FINDING ANGLE (GETS ANGLE FROM RAW XY FOR EACH LED)
%     ExtractedAngle=XYangle(data_video(:,2),data_video(:,3));
    ExtractedAngle=XYangleLED(table(:,2),table(:,4),table(:,3),table(:,5));
    
    data_video=[data_video,ExtractedAngle,[vel_cmPerSec; vel_cmPerSec(end)]]; % remove first point to fit with vel
    
    % EXTRACT BASIC MOVEMENT DATA
    BasicLoco.AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*sampleRate));
    BasicLoco.OverallDistTraveled=sum(vel_abs*pixelDist);
    BasicLoco.MeanVelocity=mean(vel_cmPerSec);
    data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).BasicLoco=BasicLoco;
    
    % START OF MAIN SPIKE TO PATH AND DATA ANALYSIS LOOP
    for i=1:length(S)
        cell=strsplit(ID{i},filesep);
        disp([cell{end},'   Cell: ',num2str(i)])
        
        SpikeFile=S{i}*1000000;
        PerSpkslessthan2ms=(sum((diff(S{i})*1000)<2)/length(S{i}))*100; %DOUBLE CHECK THIS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
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
        TS=interp1(data_video_nospk(:,1),data_video_nospk(:,1),SpikeFile,'linear');
        X=interp1(data_video_nospk(:,1),data_video_nospk(:,2),SpikeFile,'linear');
        Y=interp1(data_video_nospk(:,1),data_video_nospk(:,3),SpikeFile,'linear');
        A=interp1(data_video_nospk(:,1),data_video_nospk(:,4),SpikeFile,'linear');
        VEL=interp1(data_video_nospk(:,1),data_video_nospk(:,5),SpikeFile,'linear');
        VELidx=interp1(data_video_nospk(:,1),data_video_nospk(:,6),SpikeFile,'linear');
        
        % CONCAT AND SORT
        data_video_spk=sortrows([[TS X Y A VEL VELidx ones(size(TS,1),1)];[data_video_nospk,zeros(length(data_video_nospk),1)]],1);
        
        % VELO FILTER BASED ON INDEX CREATED ABOVE
        data_video_nospk(logical(in),:)=[];
        data_video_nospk(:,6)=[];
        nspkbefore=sum(data_video_spk(:,7)==1);
        data_video_spk(data_video_spk(:,6)==1,:)=[];
        nspkafter=sum(data_video_spk(:,7)==1);
        data_video_spk(:,6)=[];
        
        disp(['Velocity Filter: ',num2str(round(sum(in)/30)),' Seconds ', 'and ',num2str(nspkbefore-nspkafter),' out of ',num2str(length(SpikeFile)),' spikes'])
      
        % IFR
        ifr=IFR(data_video_spk(data_video_spk(:,6)==1,1),data_video_spk(:,6),round(length(data_video_nospk)/30),30);
        VIFR=[data_video_spk(:,5),ifr];
        VIFR(isnan(VIFR),:)=[];
        SpeedScore=corr(VIFR(:,1),VIFR(:,2));
        
        % THETA MODULATION
        [thetaindex,thetapeak,cor,~]=thetamodulation(data_video_spk(data_video_spk(:,6)==1,1)); %LOOK INTO HOW THIS WORKS
      
            spks_VEL=data_video_spk(data_video_spk(:,6)==1,:); %ALL FRAME DATA THAT HAS SPIKES (FOR BINNING)
            
            %%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%
            % BIN & SMOOTH CIRCLE DATA
                [SmoothRateMap,nBinsx,nBinsy,occ,Coherence]=bindata(data_video_nospk,sampleRate,spks_VEL,linear_track,track_length);
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
                else
                    infieldFR=0;outfieldFR=0;
                end
                
                % CALCULATE LANDMARK CONTROL FOR PLACE CELL
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
                    
                    % rescale xy to length of bins
                    frames=[rescale(data_video_spk(:,2),1,length(SmoothRateMap)),rescale(data_video_spk(:,3),1,length(SmoothRateMap))];
                    
                    % find durations field at least 200ms and over
                    in=inpolygon(frames(:,1),frames(:,2),bound(:,1),bound(:,2));
                    dsig=diff([0 (abs(in')>=eps) 0]);
                    startIndex=find(dsig>0);endIndex=find(dsig<0)-1;
                    stringIndex=(endIndex-startIndex+1>=6);startIndex=startIndex(stringIndex);endIndex=endIndex(stringIndex);
                    indices=zeros(1,max(endIndex)+1);indices(startIndex)=1;indices(endIndex+1)=indices(endIndex+1)-1;indices=find(cumsum(indices));
                    in=zeros(length(in),1);in(indices',1)=1;
                    pass=find(diff([0 indices])>1==1); pass=[1,pass];
                    % normalize passes from 0 to 1
                    normalizedD=zeros(length(in),1);
                    for j=1:length(pass)
                        if j+1>length(pass);break;end
                        tempidx=indices(pass(j):pass(j+1)-1);
                        normalizedD(tempidx',1)=linspace(0,1,length(tempidx))';
                    end
                    
                    % index out times in field that contain spikes
                    spks_VEL4LFP=data_video_spk(logical(in) & data_video_spk(:,6)==1,:);
                    normalizedD=normalizedD(logical(in) & data_video_spk(:,6)==1,:);
                else
                    spks_VEL4LFP=data_video_spk(data_video_spk(:,6)==1,:);
                    normalizedD=linspace(0,1,size(spks_VEL4LFP,1))';
                    E=NaN;
                end
                ~~~~~~~~~~~~~~~~~~~~
            % LFP ANALYSIS
            % Create name of CSC file from cell and path name
            %             [~, filename] = fileparts(tfile{i}); tetrode=regexp(filename,'.','match');eegfile=cellstr(strcat(path,filesep,'CSC',tetrode(3),'.ncs'));
            [~, filename] = fileparts(ID{i}); tetrode=regexp(filename,'.','match');eegfile=cellstr(strcat(path,filesep,'CSC',tetrode(end),'.ncs'));
            
            [ThPrecess,ThetaStats]=EEGWorking2(eegfile{1},spks_VEL4LFP,StartofRec,EndofRec,event,track_length,normalizedD);
            
            % spikes in field
            nspikes_infield=sum(spks_VEL4LFP(:,6));
            
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
                    fig=figure;subplot(5,1,1),plot(plotframes(:,2),plotframes(:,3),'LineWidth',1,'color','k');
                    hold on;box off; axis off
                    scatter(spks_VEL(:,2),spks_VEL(:,3),20,'filled','r');
                    title([cell{end},' Cell: ',num2str(i),'  nSpikes: ',num2str(nSpikes)]);
                    
                    % plot filled and smoothed rate map
                    imAlpha=ones(size([SmoothRateMap;SmoothRateMap]));
                    imAlpha(isnan([SmoothRateMap;SmoothRateMap]))=0;
                    subplot(5,1,2); imagesc([SmoothRateMap;SmoothRateMap],'AlphaData',imAlpha);
                    hold on; shading interp; colormap jet; axis off; box off;axis image
                    title(['InfoContent: ',num2str(round(InformationContent,3)),' F2W: ',num2str(round(Field2Wall,3))]);
                    
                    subplot(5,1,3), area(SmoothRateMap(1,:),'LineWidth',2,'EdgeColor',[0,0,0]+0.4,'FaceColor',[0,0,0]+0.8);
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
                        subplot(5,1,4); plot(ThPrecess.scatteredPH(:,1),ThPrecess.scatteredPH(:,2),'k.');hold on;
                        plot(ThPrecess.scatteredPH(:,1),ThPrecess.scatteredPH(:,2)+360,'r.')
                        ylim([0 720]); xlim([0 1])
                        set(gca, 'YTick', [0;240;480;720],'Box','off');
                        title(['Precess R^2: ',num2str(round(ThPrecess.RSquared,3)),', Slope: ',...
                            num2str(round(ThPrecess.slope,3)),', Corr: ',num2str(round(ThPrecess.Correlation,3))])
                        
                        % PLOT SMOOTHED RATEMAP
                        subplot(5,1,5); h=pcolor(ThPrecess.smoothedPHmap);
                        shading interp; colormap jet; hold on; box off; axis off; set(h, 'EdgeColor', 'none');
                        title(['DOM: ',num2str(round(ThPrecess.DOM,3))])
                        
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
            
            % pack results for NaN measures given maze type
            if isequal(linear_track,'no');DirectionalityIndex=NaN; else; borderScore=NaN; E=NaN;DisplacementCorr=NaN;end
            if ~exist('Displacement','var');Displacement=NaN;DisplacementCorr=NaN;end
            
            measures=[InformationContent,Coherence,sparsity,PeakRate,...
                OverallFR,Field2Wall,FieldWidth,nSpikes,mean_vector_length,...
                preferred_Direction(1),Direct_infoContent,DirectionalityIndex,...
                ThPrecess.Correlation,ThPrecess.DOM,ThPrecess.slope,ThPrecess.meanFR,ThPrecess.RSquared,...
                ThetaStats.MeanThetaFreq,ThetaStats.MeanOverallPow,ThetaStats.ThetaRatio,...
                PerSpkslessthan2ms,Displacement,infieldFR,outfieldFR,grade(i)...
                borderScore,E,NumbActiveBins,DisplacementCorr,SpeedScore,IntraTrialR,nspikes_infield,...
                ThPrecess.phaselock.Rlength,ThPrecess.phaselock.Pval,thetaindex,thetapeak];
            
            %SAVE SESSION DATA
            data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,event)=measures;
            data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ratemap.(['Cell',num2str(i)]).(['session',num2str(event)])=SmoothRateMap;
            data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).thetaautocorr.(['Cell',num2str(i)]).(['session',num2str(event)])=cor;
            
      
        keep('path','data_video','SpikeFile','mclustpath','tfile','track_length',...
            'linear_track','counter','SmoothRateMap_Right','SmoothRateMap_Left'...
            ,'figures','S','sampleRate','i','Overall_DirectionalityIndex','StartofRec',...
            'EndofRec','event','timestamps','vel_cmPerSec','table2',...
            'lowercut','highercut','pixelDist','grade','data','ratID','video');
    end
end

% SAVE RESULTS TO MAT FILE
if exist('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','file')>0
    datatemp=data;
    load('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','data')
    data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')])=datatemp.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]);
end
save('D:\Place_Cell_Data\RawPAE_PlaceCell\data.mat','data','-v7.3')

disp 'DONE'

%~~~~~~~~~~~~~~CODE GRAVEYARD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Arena diameter - INCLUDE LATER LB 11-Feb-2018 (see code graveyard)
%     ratID=strsplit(path,filesep);
%     if isequal(linear_track,'no') && strcmp(ratID(end-1),'LE2813') || strcmp(ratID(end-1),'LE2821') && event==4
%         track_length = 100;
%     elseif isequal(linear_track,'no')
%         track_length = 76.5;
%     end
