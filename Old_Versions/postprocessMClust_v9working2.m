
%%%%%%%% DATA EXTRACTION & POSTPROCESSING from Neuralynx Files %%%%%%%%%%%%%
%
% THIS SCRIPT CREATES FIGURES AND STATS FOR PLACE & HD CELL DATA
% CAN HANDLE CIRCULAR ARENA & LINEAR TRACK DATA
%
% INPUTS: PATH TO RAW NON-PROCESSED DATA (NON-P FOLDER) / TRACK LENGTH / MAZE-TYPE / FIGURES
%
% OUTPUTS: RATE MAP / SPIKE ON PATH / TUNING CURVE / POLAR PLOT / R VS. L SPIKES / STATS / LFP ANALYSIS
% Ryan Harvey & Ben Clark Updated: 9/21/17

%    __    __    __    __    __    __    __    __    __    __    __
% __/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__
%   \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/


function postprocessMClust_v9working2(path,track_length,linear_track,figures)
close all;clc
warning('off', 'MATLAB:dispatcher:nameConflict')
% ADD TOOLBOXES TO PATH
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));  % PC
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analyses\spikeCode'));  % PC
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\chronux_2_11'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\chronux_2_11'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\BClarkToolbox\Analysis'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\CircStat2012a'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\CircStat2012a'));
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\FMAToolbox'));
addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\FMAToolbox'));
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analyses/spikeCode'));  % MAC
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/chronux_2_11'));
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis'));
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/CircStat2012a'));         
addpath(genpath('/Users/RyanHarvey/GoogleDrive/MatlabDir/FMAToolbox'));
% addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/images'));

% SPLIT UP SESSIONS BY EVENT FILES
[ timestamps,StartofRec,EndofRec ] = EventSplit( path );                   
if isempty(timestamps)
    clear timestamps
end

% START OF EVENT LOOP 1:nSESSIONS
for event=1:length(StartofRec)
    if event==1 && isequal(linear_track,'yes'); linear_track = 'yes'; else linear_track = 'no'; end
    
    % READ VIDEO FILE
    % Creates video tracking text file
    if exist([path,filesep,'VT1.txt'],'file')==0; ReadVideoTrackingFile3(path); end
    
    % locates newly created VT1.txt file
    VTfile = FindFiles('VT1.txt', 'StartingDirectory',path);
    table=textread(cell2mat(VTfile)); % reads table
    
    % RESTRICT DATA BY START AND END OF EVENT
    if exist('timestamps','var')==1;
        try
            table=table(table(:,1)>StartofRec(event) & table(:,1)<EndofRec(event),:);
        catch
            disp(['Event File Failure in Session',num2str(event)]);
            continue
        end
    end
    
    % NEW LENGTH OF DATA            30 TIMESTAMPS = 1 SECOND
    lengthofdata=length(table);     sampleRate = 30;
    
    % DURATION OF SESSION (MIN)
    sessionduration=(lengthofdata/sampleRate)/60;
    disp(['SESSION_',num2str(event),' WAS ',num2str(sessionduration),' MIN'])
    if sessionduration<2;
        disp('SESSION TOO SHORT...SKIPPING')
        continue
    end
    
    % DEAL WITH NON DETECTS
    [table,SumofNaN1,SumofNaN2,~] = NonDetects(table,lengthofdata);
    
    % add NSMA to path for FindFiles
    addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
    addpath(genpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
    addpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox'));
    
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
    if isequal(linear_track,'no')
        table = table((vel_cmPerSec<=100),:);
    else
        table = table((vel_cmPerSec<=100),:);
    end
    
    filteredlengthofdata=length(table);
    disp(['VELOCITY FILTERING ', num2str(lengthofdata-filteredlengthofdata),...
        ' DATA POINTS, WHICH IS ', num2str(100-((filteredlengthofdata/lengthofdata)*100)), ' PERCENT OF YOUR DATA.']);
    
    % SMOOTHING RAW XY DATA FROM EACH LED
    padDsData=[repmat(table(1,:),30,1); table; repmat(table(end,:),30,1)]; %pad ends with last data point
    padDsData=[padDsData(:,1),runline(padDsData(:,2),10,1),runline(padDsData(:,3),10,1),runline(padDsData(:,4),10,1),runline(padDsData(:,5),10,1)];
    table=(padDsData(30+1:end-30,:)); %remove pad
    %     table=[table(:,1),runline(table(:,2),10,1),runline(table(:,3),10,1),runline(table(:,4),10,1),runline(table(:,5),10,1)];
    
    % GET VELOCITY
    [vel_cmPerSec,vel_abs,pixelDist] = InstaVel([table(:,VelIndex(1)),table(:,VelIndex(2))],linear_track,track_length,sampleRate);
    
    % CONCAT TIMESTAMP/XY POSITION/ANGLE TOGETHER
    data_video_smoothfilt2=double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2)]);
    
    % STRAIGHTEN LINEAR TRACK
    if isequal(linear_track,'yes')
        [X,Y] =ParameterizeToLinearTrack2(data_video_smoothfilt2(:,2),data_video_smoothfilt2(:,3));
        data_video_smoothfilt2=[data_video_smoothfilt2(:,1) X Y];
    end
    
    % FINDING ANGLE (GETS ANGLE FROM RAW XY FOR EACH LED)
    ExtractedAngle=XYangle(data_video_smoothfilt2(:,2),data_video_smoothfilt2(:,3));
    
    data_video_smoothfilt2=[data_video_smoothfilt2(2:end,:),ExtractedAngle,vel_cmPerSec]; % remove first point to fit with vel
    
    % EXTRACT BASIC MOVEMENT DATA
    BasicLoco.AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*sampleRate));
    BasicLoco.OverallDistTraveled=sum(vel_abs*pixelDist);
    BasicLoco.MeanVelocity=mean(vel_cmPerSec);
    
    % sets path to p file & looks within TT folder for tfiles
    mclustpath=cd([path,'p',filesep 'TT']);
    
    ls '*.t' % disp .t files
    tfile=FindFiles('*.t');
    
    counter = 1;
    % loads .t into ts cells
    S =LoadSpikes(tfile);
    
    % START OF MAIN SPIKE TO PATH AND DATA ANALYSIS LOOP
    for i=1:length(tfile) % Counts # of .t files (# of clusters) and sets iteration length
        warning('off', 'MATLAB:dispatcher:nameConflict')
        addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
        addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
        addpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox')); % add NSMA to path
        
        % DISPLAY FILE, CLUSTER, & ITERATION
        disp([tfile{i},'    ' num2str(counter), ' of ', num2str(length(tfile))]);
        
        % Calculate percent of spikes that occur less than 2ms & 60hz
        [filepath,filename]=fileparts(tfile{i}); [TTclust]= strsplit(filename,'_');
        if exist(char(strcat(filepath,filesep,TTclust(1), '_ClusterSummary_', TTclust(2),'.mat')),'file')>0
            load(char(strcat(filepath,filesep,TTclust(1), '_ClusterSummary_', TTclust(2),'.mat')),'CI');
            PerSpkslessthan2ms=(sum(CI.HistISI(CI.HistISICenters<=2))/CI.nSpikes)*100;
            newCorr=CI.AutoCorr-mean(CI.AutoCorr); newCorr(newCorr<0)=0; % 60hz check
        else
            PerSpkslessthan2ms=NaN;
            newCorr=NaN;
        end
        
%         if PerSpkslessthan2ms>3 || sum(newCorr>0)==60 % sum(newCorr>0)<62 && sum(newCorr>0)>60
%             [counter]=deletefiltered(tfile{i},counter,event);
%             disp(['CELL ',filename, ' ISI <2ms >3%, or 60hz detected'])
%             continue
%         end
        
        % conversion from MClust output to our current code
        SpikeFile=100*Data(S{i});
        
%         if size(SpikeFile,1)<50
%             [counter]=deletefiltered(tfile{i},counter,event);
%             disp(['CELL ',filename, ' <50 SPIKES'])
%             continue
%         end
        
        % RESTRICT DATA BY START AND END OF EVENT
        if exist('timestamps','var')==1;
            SpikeFile=SpikeFile(SpikeFile(:,1)>StartofRec(event) & SpikeFile(:,1)<EndofRec(event),:);
        end
        
        rmpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)
        rmpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)
        rmpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)
        
        % INTERPOLATE SPIKES TO TIMESTAMPS, POSITION, AND VEL DATA
        TS = interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,1), SpikeFile, 'linear');
        X= interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,2), SpikeFile, 'linear');
        Y = interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,3), SpikeFile, 'linear');
        A = interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,4), SpikeFile, 'linear'); % *** NEED TO USE DIFFERENT INTERP FOR ANGLULAR DATA ***
        VEL = interp1(data_video_smoothfilt2(:,1), data_video_smoothfilt2(:,5), SpikeFile, 'linear');
        
        % CONCAT AND SORT
        data_video_smoothfilt=sortrows([[TS X Y A VEL ones(size(TS,1),1)];[data_video_smoothfilt2,zeros(length(data_video_smoothfilt2),1)]],1);
        
        % VEL FILTER FOR MOVEMENT
%         data_video_smoothfilt3=VelFilter( data_video_smoothfilt2,data_video_smoothfilt2(:,5),2 );
%         data_video_smoothfilt=VelFilter( data_video_smoothfilt,data_video_smoothfilt(:,5),2 );
        
        % CUT OFF ENDS OF TRACK
        if strcmp(linear_track,'yes')==1
            [data_video_smoothfilt3,data_video_smoothfilt,track_length2] = RemoveEnds(pixelDist,data_video_smoothfilt2,data_video_smoothfilt,track_length,5);
        end
        
        % SPLIT UP DATA BY RIGHT AND LEFT RUNS
        if isequal(linear_track, 'yes'); runiteration=2; else runiteration=1; end
        if isequal(linear_track, 'yes')
            [right,left,DirectionalityIndex,Displacement]=RightVsLeft(data_video_smoothfilt3,data_video_smoothfilt,track_length2,sampleRate);
        end
        for iruns=1:runiteration
            % UNPACK STRUCTURE
            if isequal(linear_track, 'yes')
                if iruns==1
                    occ=right.occ; nBinsx=right.nBinsx; nBinsy=right.nBinsy; Coherence=right.Coherence;
                    SmoothRateMap=right.SmoothRateMap; data_video_smoothfilt3=right.datavid; data_video_smoothfilt=right.dataspks;
                elseif iruns==2
                    occ=left.occ; nBinsx=left.nBinsx; nBinsy=left.nBinsy; Coherence=left.Coherence;
                    SmoothRateMap=left.SmoothRateMap; data_video_smoothfilt3=left.datavid; data_video_smoothfilt=left.dataspks;
                end
                Displacement=Displacement*(track_length/length(SmoothRateMap));
            end
            spks_VEL = data_video_smoothfilt(data_video_smoothfilt(:,6) == 1,:);
            
%             % filter by less than 5 spikes
%             if size(spks_VEL,1)<50 && isequal(linear_track, 'no')
%                 clear 'spks_VEL' 'data_video_smoothfilt' % REMOVE SPIKES
%                 [counter]=deletefiltered(tfile{i},counter,event);
%                 disp(['CELL ',filename, ' LESS THAN 50 SPIKES'])
%                 continue
%             end
            
            %%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%
            
            % LFP ANALYSIS
            % Create name of CSC file from cell and path name
            [~, filename] = fileparts(tfile{i}); tetrode=regexp(filename,'.','match'); eegfile=cellstr(strcat(path,filesep,'CSC',tetrode(4),'.ncs'));
            
            % Compute Phase Lock
            if isequal(linear_track,'yes')
                % ANALYZE JUST FIRING FIELD
                [spks_VEL4LFP,FieldWidth,field]=FindFF(SmoothRateMap,data_video_smoothfilt,pixelDist,track_length2);
                
                Map4bordertemp=SmoothRateMap;Map4bordertemp(~field)=NaN;
                infieldFR=nanmean(reshape(Map4bordertemp,[],1));
                
                Map4bordertemp=SmoothRateMap;Map4bordertemp(logical(field))=NaN;
                outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
                
                if size(spks_VEL,1)>=2
                    [Stats,ThPrecess,ThetaStats]=EEGWorking2(eegfile{1},spks_VEL,SpikeFile,StartofRec,EndofRec,event,linear_track,track_length2);
                end
            else
                if size(spks_VEL,1)>=2
                    [Stats,ThPrecess,ThetaStats]=EEGWorking2(eegfile{1},spks_VEL,SpikeFile,StartofRec,EndofRec,event,linear_track,track_length);
                end
            end
            
            % BIN & SMOOTH DATA
            if isequal(linear_track,'no')
                [SmoothRateMap,nBinsx,nBinsy,occ,Coherence] = bindata(data_video_smoothfilt3,sampleRate,spks_VEL,linear_track,track_length);
            end
            
            % calculate information content
            rY = reshape(SmoothRateMap,nBinsx*nBinsy,1); % reshape data into column
            rY(isnan(rY))=0; rY(isinf(rY))=0;
            occRSHP = reshape(occ,nBinsx*nBinsy,1); % reshape data into column
            occSUM = sum(occRSHP); % summed occupancy
            pX = occRSHP./occSUM; % normalized occupancy
            [InformationContent] = InformationPerSpike(rY,pX); % from NSMA toolbox
            
            % Calculate the distance from place field to wall - peaks method
            r_max = max(rY);
            PeakRate = r_max;
            [~, x_max] = find(SmoothRateMap == r_max);
            x_top = nBinsx - x_max;
            x_bottom = x_max - 0;
            distanceAll = [x_top x_bottom];
            distanceMin = min(distanceAll);
            DistanceFromTrackEnd = distanceMin*(track_length/length(SmoothRateMap));
            
            % OVERALL FIRING RATE
            OverallFR = (size(spks_VEL,1)/size(data_video_smoothfilt3,1))*sampleRate;
            
            % NUMBER OF SPIKES
            nSpikes=size(spks_VEL,1);
            
            % calculate sparsity
            [sparsity] = Sparsity(rY'); % from NSMA toolbox
            
            % calculate percent of active bins (estimate of field size - Royer et al., 2010)
            NumbActiveBins = numel(find(rY > r_max*0.20));
            
            if isequal(linear_track,'no')
                SmoothRateMap2=SmoothRateMap;
                [field,FieldWidth]=FindFF2D(SmoothRateMap2);
                FieldWidth=FieldWidth*(track_length/length(field));
                if isnan(field)==0
                    nanPlacement=isnan(SmoothRateMap2);
                    SmoothRateMap2(~field)=0; SmoothRateMap2(nanPlacement)=NaN;
                    
                    Map4bordertemp=SmoothRateMap2;Map4bordertemp(~field)=NaN;
                    infieldFR=nanmean(reshape(Map4bordertemp,[],1));
                    
                    Map4bordertemp=SmoothRateMap;Map4bordertemp(logical(field))=NaN;
                    outfieldFR=nanmean(reshape(Map4bordertemp,[],1));
                else
                    infieldFR=0;
                    outfieldFR=0;
                end
                % find field
                [row,col]=find(field);
                %find boundary
                k=boundary(row,col);
                bound(i).control=[col(k),row(k)];
                E=Eccentricity(col(k),row(k));
                try
                    [borderScore,Field2Wall]=BorderScore(SmoothRateMap2);
                    Field2Wall=Field2Wall*(track_length/length(SmoothRateMap2));
                catch
                    borderScore=NaN;
                    Field2Wall=NaN;
                end
            end
            
            % GRID CELL ANALYSIS
            % GRID SCORE
            if isequal(linear_track,'no');
                try
                    autocorr=SpatialAutoCorr(SmoothRateMap,length(SmoothRateMap));
                    gridout=GridScore_Sinusoidal(autocorr,length(SmoothRateMap));
                    gridscore=gridout.maxSinuGrid;
                catch
                    gridscore=NaN;
                end
            end
            % SPIKE DIRECTION (from Shawn Winter's code)
            [mean_vector_length,peak_Firing_Rate,preferred_Direction,halfPeak,Directional_Range_HalfWidth,Direct_infoContent,BinsNbSpikes,BinsAngle3,BinsAngle]...
                = HDCell(data_video_smoothfilt(:,6),data_video_smoothfilt(:,4),sampleRate);
            
            % -------------------------------CREATING PLOTS----------------------------
            if figures==1
                if isequal(linear_track, 'yes')
                    % plot spike on path
                    %                     figure (i), subplot(5,1,1), plot(data_video_smoothfilt(:,2), data_video_smoothfilt(:,3), 'LineWidth', 1, 'color', 'k');
                    
                    [vel_,~,~] = InstaVel(data_video_smoothfilt,'no',track_length,sampleRate);
                    point2=find(vel_>prctile(data_video_smoothfilt(:,5),90)==1);
                    point=1;
                    figure(i); subplot(5,1,1)
                    for ii = 1:length(point2)
                        plot(data_video_smoothfilt(point:point2(ii),2),data_video_smoothfilt(point:point2(ii),3), 'LineWidth', 1, 'color', 'k')
                        hold on; point=point2(ii)+1;
                    end
                    hold on; scatter(spks_VEL(:,2), spks_VEL(:,3), 20, 'filled', 'r'); box off; axis off
                    title(['Spike on Path /',' nSpikes: ',num2str(nSpikes)]);
                    
                    % plot filled and smoothed rate map
                    figure (i), subplot(5,1,2), h = pcolor([SmoothRateMap;SmoothRateMap]);
                    hold on; shading interp; colormap jet; axis off
                    hold on; box off; set(h, 'EdgeColor', 'none'); axis image
                    title(['Smoothed RateMap',' InfoContent: ',num2str(InformationContent)]);
                    
                    
                    figure (i), subplot(5,1,3), area(SmoothRateMap(1,:),'LineWidth',2,'EdgeColor',[0,0,0]+0.4,'FaceColor', [0,0,0]+0.8);
                    box off; xlim([1 nBinsx]);
                    title(['Rate Plot',', Peak Rate: ',num2str(PeakRate)]);
                    set(figure (i),'Position',[842 42 838 924]);
                    
                    % PLOT SCATTERED PHASE PRECESSION
                    if exist('ThPrecess','var')
                        trackX=ThPrecess.scatteredPH(:,1).*pixelDist;
                        figure(i); subplot(5,1,4); plot(trackX, ThPrecess.scatteredPH(:,2), 'k.')
                        hold on; plot(trackX, ThPrecess.scatteredPH(:,2)+360, 'r.')
                        ylim([0 720]); xlim([min(data_video_smoothfilt(:,2))*pixelDist,max(data_video_smoothfilt(:,2))*pixelDist])
                        set(gca, 'YTick', [0;240;480;720],'Box','off');
                        title(['PhasePrecession',' R-Squared: ',num2str(ThPrecess.RSquared),...
                            ', Slope: ',num2str(ThPrecess.slope),', Corr: ',num2str(ThPrecess.Correlation)])
                        hold off
                        
                        % PLOT SMOOTHED RATEMAP
                        figure(i); subplot(5,1,5); h=pcolor(ThPrecess.smoothedPHmap);
                        shading interp; colormap jet; hold on; box off; axis off; set(h, 'EdgeColor', 'none');
                        title(['Mean Firing Rate: ',num2str(ThPrecess.meanFR),', DOM: ',num2str(ThPrecess.DOM)])
                    end
                    fig1 = figure(i);
                    
                    % FOR CIRCULAR ARENA
                elseif isequal(linear_track, 'no')
                    % plot spike on path
                    figure (i), subplot(2,2,1); plot(data_video_smoothfilt(:,2), data_video_smoothfilt(:,3), 'LineWidth', 1, 'color', 'k');
                    hold on; axis off
                    scatter(spks_VEL(:,2), spks_VEL(:,3), 35, 'filled', 'r');
                    box off; axis image
                    title(['Spike on Path, nSpikes: ',num2str(nSpikes)]);
                    
                    
                    figure (i), subplot(2,2,2); h = pcolor(SmoothRateMap);
                    colormap jet; axis off; hold on; box off; set(h, 'EdgeColor', 'none'); axis image; %shading interp;
                    title(['Smoothed Rate Map, IC: ',num2str(InformationContent)]);
                    
                    %smooth spike data
                    [pdf,~]=circ_ksdensity(spks_VEL(:,4), 0:359,'msni');
                    
                    % Plot Smoothed firing rate x HEAD DIRECTION
                    figure (i),subplot(2,2,3); plot(pdf,'LineWidth',2,'color','k')
                    axis tight; hold on; xlim ([0 360]); box off
                    title(['Tuning Curve, Peak Rate: ',num2str(peak_Firing_Rate),' D_IC: ',num2str(Direct_infoContent)])
                    
                    % Firing rate x HD polar plot for the nonsmoothed data above
                    figure(i); subplot(2,2,4); polarplot = polar(BinsAngle([1:60 1]),BinsNbSpikes([1:60 1]),'b');
                    set(polarplot, 'linewidth',3,'color','k'); axis off
                    title(['Polor Plot, MeanVecLength: ',num2str(mean_vector_length),' Pref_Dir: ',num2str(preferred_Direction)]);
                    set(0,'Showhiddenhandles','on')
                    set(figure(i),'Position',[686 325 977 619]);
                    fig1 = figure(i);
                end
            end
            
            % SAVING JPGS AND .MAT FILES TO TT FOLDER
            cd(mclustpath);
            warning('off', 'MATLAB:Figure:FigureSavedToMATFile');
            % SAVE RUN ITERATION
            if runiteration==2
                % DELETE OLD MAT FILE FROM BEFORE SPLITING RUNS
                [~]=deletefiltered(tfile{i},counter,event);
                % SAVE WITH S# IF EVENTS EXIST
                if exist('timestamps','var')
                    [filepath, filename] = fileparts(tfile{i});
                    if figures==1
                        saveas(fig1,[filepath filesep filename sprintf('R%d',iruns) sprintf('S%d',event) '_spikeFigure.jpg']);
                    end
                    save([filepath filesep filename sprintf('R%d',iruns) sprintf('S%d',event) '_spikeData.mat']);
                else
                    [filepath, filename] = fileparts(tfile{i});
                    if figures==1
                        saveas(fig1,[filepath filesep filename sprintf('R%d',iruns) '_spikeFigure.jpg']);
                    end
                    save([filepath filesep filename sprintf('R%d',iruns) '_spikeData.mat']);
                end
                close all
            else
                % SAVE WITH S# IF EVENTS EXIST
                if exist('timestamps','var');
                    [filepath, filename] = fileparts(tfile{i});
                    if figures==1
                        saveas(fig1,[filepath filesep filename sprintf('S%d',event) '_spikeFigure.jpg']);
                    end
                    save([filepath filesep filename sprintf('S%d',event) '_spikeData.mat']);
                else
                    [filepath, filename] = fileparts(tfile{i});
                    if figures==1
                        saveas(fig1,[filepath filesep filename '_spikeFigure.jpg']);
                    end
                    save([filepath filesep filename '_spikeData.mat']);
                end
                close all
            end
        end
        keep('path','data_video_smoothfilt2','SpikeFile','mclustpath','tfile','track_length',...
            'linear_track','counter','SmoothRateMap_Right','SmoothRateMap_Left'...
            ,'figures','S','sampleRate','i','Overall_DirectionalityIndex','StartofRec',...
            'EndofRec','event','timestamps','vel_cmPerSec','table2',...
            'lowercut','highercut','pixelDist');
        counter = counter + 1;
    end
end
disp 'DONE'




