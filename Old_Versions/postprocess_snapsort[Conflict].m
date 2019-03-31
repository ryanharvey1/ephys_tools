
%%%%%%%% DATA EXTRACTION & POSTPROCESSING from Neuralynx SNAP sort file %%%%%%%%%%%%%
%
% THIS SCRIPT CREATES FIGURES AND STATS FOR PLACE & HD CELL DATA
% CAN HANDLE CIRCULAR ARENA & LINEAR TRACK DATA
%
% INPUTS: PATH TO RAW NON-PROCESSED DATA (NON-P FOLDER) / TRACK LENGTH / MAZE-TYPE / FIGURES
%
% OUTPUTS: RATE MAP / SPIKE ON PATH / TUNING CURVE / POLAR PLOT / R VS. L SPIKES / STATS / LFP ANALYSIS
% Ryan Harvey: 12/14/17

%    __    __    __    __    __    __    __    __    __    __    __
% __/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__
%   \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/


function postprocess_snapsort(path,track_length,linear_track,figures)
close all;clc
% ADD TOOLBOXES TO PATH
cd F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox
com=pwd;com=strsplit(com,filesep);com=com{1};
addpath([com,'\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode\MEX']);  % PC
addpath([com,'\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode']);  % PC
addpath([com,'\Users\BClarkLab\GoogleDrive\MatlabDir\chronux_2_11']);
addpath([com,'\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis']);
addpath([com,'\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis\LFP'])
addpath([com,'\Users\BClarkLab\GoogleDrive\MatlabDir\CircStat2012a']);
addpath([com,'\Users\BClarkLab\GoogleDrive\MatlabDir\FMAToolbox']);

% SPLIT UP SESSIONS BY EVENT FILES
[ timestamps,StartofRec,EndofRec ]=EventSplit(path);if isempty(timestamps);clear timestamps;end

% SET UP VARIABLE NAMES & DATA STRUCTURE
varnames={'InformationContent','Coherence','Sparsity','PeakRate',...
    'OverallFiringRate','DistanceFromTrackEnd','FieldWidth','nSpikes','mean_vector_length',...
    'preferred_Direction','Direct_infoContent','DirectionalityIndex',...
    'RLength_Delta','RLength_Theta','RLength_Alpha','RLength_Beta','RLength_Gamma','RLength_HighGamma',...
    'Rayleigh_Delta','Rayleigh_Theta','Rayleigh_Alpha','Rayleigh_Beta','Rayleigh_Gamma','Rayleigh_HighGamma',...
    'RayleighZ_Delta','RayleighZ_Theta','RayleighZ_Alpha','RayleighZ_Beta','RayleighZ_Gamma','RayleighZ_HighGamma',...
    'MeanPhaseDeg_Delta','MeanPhaseDeg_Theta','MeanPhaseDeg_Alpha','MeanPhaseDeg_Beta','MeanPhaseDeg_Gamma','MeanPhaseDeg_HighGamma',...
    'median_Delta','median_Theta','median_Alpha','median_Beta','median_Gamma','median_HighGamma'...
    'PhPrecessCorr','DepthofModulation','PhPrecessSlope','PhPrecessMeanFR','PH_RSquared',...
    'MeanThetaFreq','MeanOverallPow','ThetaRatio','PerSpkslessthan2ms','Displacement',...
    'infieldFR','outfieldFR'};
ratID=strsplit(path,filesep);
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).varnames=varnames;


% START OF EVENT LOOP 1:nSESSIONS
for event=1:length(StartofRec)
    Session=event;
    if event==1 && isequal(linear_track,'yes');linear_track='yes';else;linear_track = 'no';end
    
    % maze type
    if  isequal(linear_track,'yes')
        mazetype='linear track';
    else
        mazetype='circular';
    end
    data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).mazetype=mazetype;

    % READ VIDEO FILE
    % Creates video tracking text file if one doesn't exist
    if exist([path,filesep,'VT1.txt'],'file')==0; ReadVideoTrackingFile3(path); end
    
    % locates newly created VT1.txt file
    fileID=fopen([path,filesep,'VT1.txt'],'r');
    dataArray=textscan(fileID,'%f%f%f%f%f%[^\n\r]','Delimiter','\t','TextType','string','EmptyValue',NaN,'ReturnOnError',false);fclose(fileID);
    table=[dataArray{1:end-1}];
    
    % RESTRICT DATA BY START AND END OF EVENT
    if exist('timestamps','var')==1;
        try
            table=table(table(:,1)>StartofRec(event) & table(:,1)<EndofRec(event),:);
        catch
            disp(['Event File Failure in Session',num2str(event)]);
            continue
        end
    end
    
    % NEW LENGTH OF DATA
    lengthofdata=length(table);sampleRate = 30;
    
    % DURATION OF SESSION (MIN)
    sessionduration=(lengthofdata/sampleRate)/60;
    disp(['SESSION_',num2str(event),' WAS ',num2str(sessionduration),' MIN'])
    if sessionduration<2;disp('SESSION TOO SHORT...SKIPPING');continue;end
    
    % DEAL WITH NON DETECTS
    [table,SumofNaN1,SumofNaN2,~]=NonDetects(table,lengthofdata);
    
    % VELOCITY FILTER & SMOOTHING
    % Arena diameter
    ratID=strsplit(path,filesep);
    if isequal(linear_track,'no') && strcmp(ratID(end-1),'LE2813') && event==4
        track_length = 100;
    elseif isequal(linear_track,'no') && strcmp(ratID(end-1),'LE2821') && event==4
        track_length = 100;
    elseif isequal(linear_track,'no')
        track_length = 76.5;
    end
    
    % CHOOSE LED WITH BEST TRACKING TO VELOCITY FILTER BY
    if SumofNaN1>SumofNaN2; VelIndex=[4 5]; else VelIndex=[2 3]; end;
    
    % CALCULATE VEL TO REMOVE JUMPS IN DATA
    [vel_cmPerSec,~,~] = InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,sampleRate);
    
    % find xy coords and spikes when rat is within velocity threshold ( <100)
    table = table((vel_cmPerSec<=100),:);
    
    filteredlengthofdata=length(table);
    disp(['VELOCITY FILTERING ', num2str(lengthofdata-filteredlengthofdata),...
        ' DATA POINTS, WHICH IS ', num2str(100-((filteredlengthofdata/lengthofdata)*100)), ' PERCENT OF YOUR DATA.']);
    
    % SMOOTHING RAW XY DATA FROM EACH LED
    padDsData=[repmat(table(1,:),30,1); table; repmat(table(end,:),30,1)]; %pad ends with last data point
    padDsData=[padDsData(:,1),smooth(padDsData(:,2),10),smooth(padDsData(:,3),10),smooth(padDsData(:,4),10),smooth(padDsData(:,5),10)];
    table=(padDsData(30+1:end-30,:)); %remove pad
    
    % GET VELOCITY
    [vel_cmPerSec,vel_abs,pixelDist]=InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,sampleRate);
    
    % CONCAT TIMESTAMP/XY POSITION/ANGLE TOGETHER
    data_video_smoothfilt2=double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2)]);
    
    % STRAIGHTEN LINEAR TRACK
    if isequal(linear_track,'yes')
        [X,Y]=ParameterizeToLinearTrack2(data_video_smoothfilt2(:,2),data_video_smoothfilt2(:,3));data_video_smoothfilt2=[data_video_smoothfilt2(:,1) X Y];
    end
    
    % FINDING ANGLE (GETS ANGLE FROM RAW XY FOR EACH LED)
    ExtractedAngle=XYangle(data_video_smoothfilt2(:,2),data_video_smoothfilt2(:,3));
    data_video_smoothfilt2=[data_video_smoothfilt2(1:end-1,:),ExtractedAngle,vel_cmPerSec]; % remove first point to fit with vel
    
    % EXTRACT BASIC MOVEMENT DATA
    BasicLoco.AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*sampleRate));
    BasicLoco.OverallDistTraveled=sum(vel_abs*pixelDist);
    BasicLoco.MeanVelocity=mean(vel_cmPerSec);
    data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).BasicLoco=BasicLoco;
    
    % sets path to snap folder & looks within for mat files
    cd([path,filesep,'SNAPSorterResults']);
    file=struct2table(dir( '**/*.mat')); t=table2cell(file(:,1));
    file=t(find(cellfun(@isempty,strfind(t,'_info.mat'))),1);
    tfile_info=t(~cellfun(@isempty,strfind(t,'_info.mat')),1);
    ns=1;grade=[];
    for m=1:length(file)
        load(file{m});load(tfile_info{m})
        grade=[grade;final_grades'];
        numspk=unique(output(:,2));
        for n=1:length(numspk);tfile{ns,1}=orig_filename;S{ns,1}=output(output(:,2)==n,1); ns=ns+1;end
    end
    disp(['    ',num2str(length(S)),' Cells'])
    counter = 1;
    
    % START OF MAIN SPIKE TO PATH AND DATA ANALYSIS LOOP
    for i=1:length(S) 
        SpikeFile=S{i}*1000000;
        PerSpkslessthan2ms=(sum((diff(S{i})*1000)<2)/length(S{i}))*100;
        
        % RESTRICT DATA BY START AND END OF EVENT
        if exist('timestamps','var')==1
            SpikeFile=SpikeFile(SpikeFile(:,1)>StartofRec(event) & SpikeFile(:,1)<EndofRec(event),:);
        end
 
        % INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
        TS=interp1(data_video_smoothfilt2(:,1),data_video_smoothfilt2(:,1),SpikeFile,'linear');
        X=interp1(data_video_smoothfilt2(:,1),data_video_smoothfilt2(:,2),SpikeFile,'linear');
        Y=interp1(data_video_smoothfilt2(:,1),data_video_smoothfilt2(:,3),SpikeFile,'linear');
        A=interp1(data_video_smoothfilt2(:,1),data_video_smoothfilt2(:,4),SpikeFile,'linear');
        VEL=interp1(data_video_smoothfilt2(:,1),data_video_smoothfilt2(:,5),SpikeFile,'linear');
        
        % CONCAT AND SORT
        data_video_smoothfilt=sortrows([[TS X Y A VEL ones(size(TS,1),1)];[data_video_smoothfilt2,zeros(length(data_video_smoothfilt2),1)]],1);
        
        % SPLIT UP DATA BY RIGHT AND LEFT RUNS FOR LINEAR TRACK DATA
        if isequal(linear_track,'yes')
            runiteration=2;
            [right,left,DirectionalityIndex,Displacement]=RightVsLeft(data_video_smoothfilt2,data_video_smoothfilt,track_length,sampleRate);
        else
            runiteration=1;
        end
        measures=[];
        for iruns=1:runiteration
            if isequal(linear_track, 'yes') % UNPACK STRUCTURE
                if iruns==1
                    occ=right.occ; nBinsx=right.nBinsx; nBinsy=right.nBinsy;
                    Coherence=right.Coherence;SmoothRateMap=right.SmoothRateMap; 
                    data_video_smoothfilttemp=right.dataspks;
                elseif iruns==2
                    occ=left.occ; nBinsx=left.nBinsx; nBinsy=left.nBinsy; 
                    Coherence=left.Coherence;SmoothRateMap=left.SmoothRateMap;
                    data_video_smoothfilttemp=left.dataspks;
                end
                [C,ia]=setdiff(data_video_smoothfilt(:,1),data_video_smoothfilttemp(:,1));
                plotframes=data_video_smoothfilt;
                plotframes(ia,:)=NaN;
                
                data_video_smoothfilt=data_video_smoothfilttemp;
                
                Displacement=Displacement*(track_length/length(SmoothRateMap));
            end
            spks_VEL=data_video_smoothfilt(data_video_smoothfilt(:,6)==1,:);
            
            %%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%
            % BIN & SMOOTH CIRCLE DATA
            if isequal(linear_track,'no')
                [SmoothRateMap,nBinsx,nBinsy,occ,Coherence]=bindata(data_video_smoothfilt2,sampleRate,spks_VEL,linear_track,track_length);
                SmoothRateMap2=SmoothRateMap;[field,FieldWidth]=FindFF2D(SmoothRateMap2);FieldWidth=FieldWidth*(track_length/length(field));
                if isnan(field)==0
                    nanPlacement=isnan(SmoothRateMap2);SmoothRateMap2(~field)=0;SmoothRateMap2(nanPlacement)=NaN;
                    Map4bordertemp=SmoothRateMap2;Map4bordertemp(~field)=NaN;infieldFR=nanmean(reshape(Map4bordertemp,[],1));
                    Map4bordertemp=SmoothRateMap;Map4bordertemp(logical(field))=NaN;outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
                else
                    infieldFR=0;outfieldFR=0;
                end
                % find field
                [row,col]=find(field);k=boundary(row,col);bound=[col(k),row(k)];
                
                % rescale xy to length of bins
                frames=[rescale(data_video_smoothfilt(:,2),1,length(SmoothRateMap)),rescale(data_video_smoothfilt(:,3),1,length(SmoothRateMap))];

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
                spks_VEL4LFP=data_video_smoothfilt(logical(in) & data_video_smoothfilt(:,6)==1,:);
                normalizedD=normalizedD(logical(in) & data_video_smoothfilt(:,6)==1,:);
                
                E=Eccentricity(col(k),row(k));
                [borderScore,Field2Wall]=BorderScore(SmoothRateMap2);
                Field2Wall=Field2Wall*(track_length/length(SmoothRateMap2));
            else
                % LINEAR
                [spks_VEL4LFP,FieldWidth,field]=FindFF(SmoothRateMap,data_video_smoothfilt,pixelDist,track_length);
                normalizedD=rescale(spks_VEL4LFP(:,2),0,1);
                Map4bordertemp=SmoothRateMap;Map4bordertemp(~field)=NaN;
                infieldFR=nanmean(reshape(Map4bordertemp,[],1));
                Map4bordertemp=SmoothRateMap;Map4bordertemp(logical(field))=NaN;
                outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
            end
            
            % LFP ANALYSIS
            % Create name of CSC file from cell and path name
            [~, filename] = fileparts(tfile{i}); tetrode=regexp(filename,'.','match');eegfile=cellstr(strcat(path,filesep,'CSC',tetrode(3),'.ncs'));
            
            [Stats,ThPrecess,ThetaStats]=EEGWorking2(eegfile{1},spks_VEL4LFP,SpikeFile,StartofRec,EndofRec,event,linear_track,track_length,normalizedD);

            % calculate information content from nsma
            rY=reshape(SmoothRateMap,nBinsx*nBinsy,1);rY(isnan(rY))=0;rY(isinf(rY))=0;occRSHP=reshape(occ,nBinsx*nBinsy,1);occSUM=sum(occRSHP);pX=occRSHP./occSUM;
            [nBins,nCells]=size(rY);relR=rY./kron(ones(nBins,1),pX'*rY);log_relR=log2(relR);log_relR(isinf(log_relR))=0;
            InformationContent=sum(kron(pX,ones(1,nCells)).*relR.*log_relR);
                        
            % Calculate the distance from place field to wall - peaks method
            r_max=max(rY);PeakRate=r_max;[~, x_max]=find(SmoothRateMap==r_max);
            DistanceFromTrackEnd=min([nBinsx-x_max x_max-0])*(track_length/length(SmoothRateMap));
            
            % OVERALL FIRING RATE
            OverallFR=(size(spks_VEL,1)/size(data_video_smoothfilt2,1))*sampleRate;
            
            % NUMBER OF SPIKES
            nSpikes=size(spks_VEL,1);
            
            % calculate sparsity from NSMA toolbox
            R=rY';if ~isa(R,'cell');R={R};end;RC=cat(2,R{:});[~,nCells]=size(RC);
            sparsity=sum(RC,2).^2./(nCells*sum(RC.^2,2));
             
            % calculate percent of active bins (estimate of field size - Royer et al., 2010)
            NumbActiveBins=numel(find(rY > r_max*0.20));
            
            % GRID CELL ANALYSIS
            % GRID SCORE
            if isequal(linear_track,'no');
                autocorr=SpatialAutoCorr(SmoothRateMap,length(SmoothRateMap));
                gridout=GridScore_Sinusoidal(autocorr,length(SmoothRateMap));
                gridscore=gridout.maxSinuGrid;
            end
            % SPIKE DIRECTION (from Shawn Winter's code)
            [mean_vector_length,peak_Firing_Rate,preferred_Direction,halfPeak,Directional_Range_HalfWidth,Direct_infoContent,BinsNbSpikes,BinsAngle3,BinsAngle]...
                = HDCell(data_video_smoothfilt(:,6),data_video_smoothfilt(:,4),sampleRate);
            
            % -------------------------------CREATING PLOTS----------------------------
            if figures==1
                if isequal(linear_track, 'yes')
                    % plot spike on path
                    figure(Session), subplot(5,1,1), plot(plotframes(:,2), plotframes(:,3), 'LineWidth', 1, 'color', 'k');
 
                    hold on; scatter(spks_VEL(:,2), spks_VEL(:,3), 20, 'filled', 'r'); box off; axis off
                    title(['Spike on Path /',' nSpikes: ',num2str(nSpikes)]);
                    
                    % plot filled and smoothed rate map
                    imAlpha=ones(size([SmoothRateMap;SmoothRateMap]));
                    imAlpha(isnan([SmoothRateMap;SmoothRateMap]))=0;
                    figure (Session);subplot(5,1,2); h=imagesc([SmoothRateMap;SmoothRateMap],'AlphaData',imAlpha);
                    hold on; shading interp; colormap jet; axis off; box off;axis image
                    title(['Smoothed RateMap',' InfoContent: ',num2str(InformationContent)]);
                    
                    figure (Session), subplot(5,1,3), area(SmoothRateMap(1,:),'LineWidth',2,'EdgeColor',[0,0,0]+0.4,'FaceColor', [0,0,0]+0.8);box off; xlim([1 nBinsx]);
                    title(['Rate Plot',', Peak Rate: ',num2str(PeakRate)]);
                    set(figure (Session),'Position',[842 42 838 924]);
                    
                    % PLOT SCATTERED PHASE PRECESSION
                    if ~isnan(ThPrecess.scatteredPH)
                        figure(Session); subplot(5,1,4); plot(ThPrecess.scatteredPH(:,1),ThPrecess.scatteredPH(:,2), 'k.') ; hold on; plot(ThPrecess.scatteredPH(:,1), ThPrecess.scatteredPH(:,2)+360, 'r.')
                        ylim([0 720]); xlim([0 1])
                        set(gca, 'YTick', [0;240;480;720],'Box','off');
                        title(['PhasePrecession',' R-Squared: ',num2str(ThPrecess.RSquared),', Slope: ',num2str(ThPrecess.slope),', Corr: ',num2str(ThPrecess.Correlation)])
                        hold off
                        
                        % PLOT SMOOTHED RATEMAP
                        figure(Session); subplot(5,1,5); h=pcolor(ThPrecess.smoothedPHmap);
                        shading interp; colormap jet; hold on; box off; axis off; set(h, 'EdgeColor', 'none');
                        title(['Mean Firing Rate: ',num2str(ThPrecess.meanFR),', DOM: ',num2str(ThPrecess.DOM)])
                    end
                    fig1 = figure(Session);
                    
                    % FOR CIRCULAR ARENA
                elseif isequal(linear_track, 'no')
                    % plot spike on path
                    figure (Session), subplot(2,2,1); plot(data_video_smoothfilt(:,2), data_video_smoothfilt(:,3), 'LineWidth', 1, 'color', 'k');
                    hold on; axis off
                    scatter(spks_VEL(:,2), spks_VEL(:,3), 35, 'filled', 'r');
                    box off; axis image
                    title(['Spike on Path, nSpikes: ',num2str(nSpikes)]);
                    
                    imAlpha=ones(size(SmoothRateMap));
                    imAlpha(isnan(SmoothRateMap))=0;
                    figure (Session);subplot(2,2,2);h=imagesc(SmoothRateMap,'AlphaData',imAlpha);
                    axis xy
                    colormap jet; axis off; hold on; box off;
                    axis image; %shading interp;
                    title(['Smoothed Rate Map, IC: ',num2str(InformationContent)]);
                    
                    %smooth spike data
                    [pdf,~]=circ_ksdensity(spks_VEL(:,4), 0:359,'msni');
                    
                    % Plot Smoothed firing rate x HEAD DIRECTION
                    figure (Session),subplot(2,2,3); plot(pdf,'LineWidth',2,'color','k')
                    axis tight; hold on; xlim ([0 360]); box off
                    title(['Tuning Curve, Peak Rate: ',num2str(peak_Firing_Rate),' D_IC: ',num2str(Direct_infoContent)])
                    
                    % Firing rate x HD polar plot for the nonsmoothed data above
                    figure(Session); subplot(2,2,4); polarplot = polar(BinsAngle([1:60 1]),BinsNbSpikes([1:60 1]),'b');
                    set(polarplot, 'linewidth',3,'color','k'); axis off
                    title(['Polor Plot, MeanVecLength: ',num2str(mean_vector_length),' Pref_Dir: ',num2str(preferred_Direction)]);
                    set(0,'Showhiddenhandles','on')
                    set(figure(Session),'Position',[686 325 977 619]);
                    fig1 = figure(Session);
                end
            end
            
            % pack results for NaN measures given maze type
            if isequal(linear_track,'no');DirectionalityIndex=NaN;Displacement=NaN;end
            
            measures=[measures;[InformationContent,Coherence,sparsity,PeakRate,...
                OverallFR,max(DistanceFromTrackEnd),FieldWidth,nSpikes,mean_vector_length,...
                preferred_Direction(1),Direct_infoContent,DirectionalityIndex,...
                Stats.MeanResultantPhase(1),Stats.MeanResultantPhase(2),Stats.MeanResultantPhase(3),...
                Stats.MeanResultantPhase(4),Stats.MeanResultantPhase(5),Stats.MeanResultantPhase(6),...
                Stats.RayleighsTest.PVal(1),Stats.RayleighsTest.PVal(2),Stats.RayleighsTest.PVal(3),...
                Stats.RayleighsTest.PVal(4),Stats.RayleighsTest.PVal(5),Stats.RayleighsTest.PVal(6),...
                Stats.RayleighsTest.ZValue(1),Stats.RayleighsTest.ZValue(2),Stats.RayleighsTest.ZValue(3),...
                Stats.RayleighsTest.ZValue(4),Stats.RayleighsTest.ZValue(5),Stats.RayleighsTest.ZValue(6),...
                Stats.Meanz(1),Stats.Meanz(2),Stats.Meanz(3),...
                Stats.Meanz(4),Stats.Meanz(5),Stats.Meanz(6),...
                Stats.median(1),Stats.median(2),Stats.median(3),...
                Stats.median(4),Stats.median(5),Stats.median(6),...
                ThPrecess.Correlation,ThPrecess.DOM,ThPrecess.slope,ThPrecess.meanFR,ThPrecess.RSquared,...
                ThetaStats.MeanThetaFreq,ThetaStats.MeanOverallPow,ThetaStats.ThetaRatio,...
                PerSpkslessthan2ms,Displacement,infieldFR,outfieldFR]];
            
            data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).ratemap.(['S',num2str(Session)])=SmoothRateMap;
            data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).frames.(['S',num2str(Session)])=data_video_smoothfilt;
            data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).Clustergrade(Session)=grade(i);
        end
        
            if isequal(linear_track,'yes') && event==1
                size(measures)
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,1)=measures(1,:);  
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,2)=measures(2,:);  
            elseif isequal(linear_track,'no') && event==1
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,1)=measures(1,:);  
            elseif event>1
                data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,event+1)=measures(1,:);  
            end
            
            
            Session=Session+1;
            
        keep('path','data_video_smoothfilt2','SpikeFile','mclustpath','tfile','track_length',...
            'linear_track','counter','SmoothRateMap_Right','SmoothRateMap_Left'...
            ,'figures','S','sampleRate','i','Overall_DirectionalityIndex','StartofRec',...
            'EndofRec','event','timestamps','vel_cmPerSec','table2',...
            'lowercut','highercut','pixelDist','grade','data','Session','measures');
        counter = counter + 1;
    end
 

                
                
%                             if isequal(linear_track,'yes') && iruns==1
%                 data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures(i,:,1)=measures;
%                 cat(3,data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).measures,measures);
%             end

end
disp 'DONE'

