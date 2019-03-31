
%%%%%%%% DATA EXTRACTION & POSTPROCESSING from Neuralynx Files %%%%%%%%%%%%%
%
% THIS SCRIPT CREATES FIGURES AND STATS FOR PLACE & HD CELL DATA
% CAN HANDLE CIRCULAR ARENA & LINEAR TRACK DATA
%
% INPUTS: PATH TO RAW NON-PROCESSED DATA (NON-P FOLDER) / TRACK LENGTH / MAZE-TYPE / FIGURES
%
% OUTPUTS: RATE MAP / SPIKE ON PATH / TUNING CURVE / POLAR PLOT / R VS. L SPIKES / STATS / LFP ANALYSIS
% Ryan Harvey & Ben Clark Updated: 2/18/17

%    __    __    __    __    __    __    __    __    __    __    __
% __/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__
%   \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/  \__/


function postprocessMClust_v9(path,track_length,linear_track,figures)
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
addpath(genpath('/Users/ryanharvey/GoogleDrive/MatlabDir/images'));

% SPLIT UP SESSIONS BY EVENT FILES
[ timestamps,StartofRec,EndofRec ] = EventSplit( path );
if isempty(timestamps)
    clear timestamps
end

% START OF EVENT LOOP 1:nSESSIONS
for event=1:length(StartofRec)
    if event==1 && isequal(linear_track,'yes'); linear_track = 'yes'; else linear_track = 'no'; end
    
    %     TimeFreqPlots(path,StartofRec,EndofRec,event)
    
    % READ VIDEO FILE
    % Creates video tracking text file
    if exist([path,filesep,'VT1.txt'],'file')==0; ReadVideoTrackingFile3(path); end
    
    % locates newly created VT1.txt file
    VTfile = FindFiles('VT1.txt', 'StartingDirectory',path);
    table=textread(cell2mat(VTfile)); % reads table
    
    % RESTRICT DATA BY START AND END OF EVENT
    if exist('timestamps','var')==1;
        table=table(table(:,1)>StartofRec(event) & table(:,1)<EndofRec(event),:);
    end
    
    % NEW LENGTH OF DATA            30 TIMESTAMPS = 1 SECOND
    lengthofdata=length(table);     sampleRate = 30; 
    
    % DURATION OF SESSION (MIN)
    sessionduration=(lengthofdata/sampleRate)/60;
    disp(['SESSION_',num2str(event),' WAS ',num2str(sessionduration),' MIN'])
    
    % DEAL WITH NON DETECTS
    [table,SumofNaN1,SumofNaN2,PercentofNaN] = NonDetects(table,lengthofdata);
    
    % VELOCITY FILTER & SMOOTHING
    % Arena diameter
    if isequal(linear_track,'no'); track_length = 76.5; end
    
    % CHOOSE LED WITH BEST TRACKING TO VELOCITY FILTER BY
    if SumofNaN1>SumofNaN2; VelIndex=[4 5]; else VelIndex=[2 3]; end;
    
    % velocity of rat from smoothed xy data
    % scalar length of velocity vector = "scalar velocity" in pixels/frame
    if isequal(linear_track,'yes');
        vel_abs = abs(diff(table(:,VelIndex(1))));
    else
        vel_abs = sqrt(diff(table(:,VelIndex(1))).^2 + diff(table(:,VelIndex(2))).^2);
    end
    pixelDist = track_length/(range(table(:,VelIndex(1))));
    vel_cmPerSec = vel_abs * pixelDist * sampleRate;
    
    % find xy coords and spikes when rat is above velocity threshold (>2 & <100)
    % Velocity filter by below 2 & above 100 cm/sec
    table2=table; % for later filtering spikes
    if isequal(linear_track,'no')
        table = table((vel_cmPerSec >=2 & vel_cmPerSec <=200),:);
    else
        table = table((vel_cmPerSec >=5),:);
    end

    filteredlengthofdata=length(table);
    disp(['VELOCITY FILTERING ', num2str(lengthofdata-filteredlengthofdata),...
        ' DATA POINTS, WHICH IS ', num2str(100-((filteredlengthofdata/lengthofdata)*100)), ' PERCENT OF YOUR DATA.']);
    
    % SMOOTHING RAW XY DATA FROM EACH LED
    table=[table(:,1),runline(table(:,2),10,1),runline(table(:,3),10,1),runline(table(:,4),10,1),runline(table(:,5),10,1)];
    
    % FINDING ANGLE (GETS ANGLE FROM RAW XY FOR EACH LED)
    ExtractedAngle = atan2d(table(:,5)-table(:,3),table(:,4)-table(:,2)) + 360*((table(:,5)-table(:,3))<0);
    
    % CONCAT TIMESTAMP/XY POSITION/ANGLE TOGETHER
    data_video_smoothfilt=double([table(:,1),median([table(:,2),table(:,4)],2),median([table(:,3),table(:,5)],2),ExtractedAngle]);
    
    % EXTRACT BASIC MOVEMENT DATA
    BasicLoco.AverageAnglePerSec=rad2deg(circ_mean((abs(diff(deg2rad(ExtractedAngle))))*sampleRate));
    BasicLoco.OverallDistTraveled=sum(vel_abs*pixelDist);
    BasicLoco.MeanVelocity=mean(vel_cmPerSec);
    
    % CUT OFF ENDS OF TRACK
    if strcmp(linear_track,'yes')==1
        data_video_smoothfilt2=data_video_smoothfilt;
        pix2remove=8/pixelDist; % Pix2remove=how many pixels in 8cms of track
        x_min=min(data_video_smoothfilt(:,2)); x_max=max(data_video_smoothfilt(:,2));
        lowercut=x_min+pix2remove; highercut=x_max-pix2remove;
        data_video_smoothfilt=data_video_smoothfilt((data_video_smoothfilt(:,2)>lowercut & data_video_smoothfilt(:,2)<highercut),:);
        track_length=track_length-16;
%         pixelDist = track_length/(range(data_video_smoothfilt(:,2)));
    end
    
    % sets path to p file & looks within TT folder for tfiles
    mclustpath=cd([path,'p',filesep 'TT']);
    
    % add NSMA to path for FindFiles
    addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
    addpath(genpath('D:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
    addpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox'));
    
    ls '*.t' % disp .t files
    tfile=FindFiles('*.t');
    
    counter = 1;
    % loads .t into ts cells
    S =LoadSpikes(tfile);
    
    % CREATE RASTER PLOT OF ALL SPIKES
    if figures==1
        RasterPlot2(S);
        [filepath, ~] = fileparts(tfile{1});
        saveas(figure(1),[filepath filesep sprintf('S%d',event) '_rasterplot.jpg']);
        close all
    end
    
    
      
    % START OF MAIN SPIKE TO PATH AND DATA ANALYSIS LOOP
    for i=1:length(tfile); % Counts # of .t files (# of clusters) and sets iteration length
        addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
        addpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
        addpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox')); % add NSMA to path
        
        % DISPLAY FILE, CLUSTER, & ITERATION
        disp([tfile{i},'    ' num2str(counter), ' of ', num2str(length(tfile))]);
        
        % conversion from MClust output to our current code
        SpikeFile=100*Data(S{i});
        
        % RESTRICT DATA BY START AND END OF EVENT
        if exist('timestamps','var')==1;
            SpikeFile=SpikeFile(SpikeFile(:,1)>StartofRec(event) & SpikeFile(:,1)<EndofRec(event),:);
        end
        
        rmpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)
        rmpath(genpath('d:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)
        rmpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)
        
        % INTERPOLATE SPIKES TO RAW PATH
        ts_spike = interp1(table2(:,1), table2(:,1), SpikeFile, 'nearest');
        
        % FILTER RAW TABLE BY VELOCITY LIMITS
        if isequal(linear_track,'no')
            table3 = table2((vel_cmPerSec >=2 & vel_cmPerSec <=200),:);
        else
            table3 = table2((vel_cmPerSec >=5),:);
        end
        
        %Create a vector of 0's and 1's (1 = ts in which spike occured) and add to array of data
        [~,ia,~] = intersect(table3(:,1), ts_spike);
        N = size(table3(:,1));
        spikeVec = zeros(N(1,1),1);
        spikeVec(ia) = 1;
        if strcmp(linear_track,'yes')==1
            spikeVec=spikeVec((data_video_smoothfilt2(:,2)>lowercut & data_video_smoothfilt2(:,2)<highercut),:);
        end
        % RECOMBINE FILTERED SPIKES WITH FILTERED AND SMOOTHED POSITION DATA
        data_video_smoothfilt = [data_video_smoothfilt spikeVec(:,1)];
                
        % Creating Plots for Right vs. Left trajectories
        if isequal(linear_track, 'yes')
            spks_VEL = data_video_smoothfilt(data_video_smoothfilt(:,5) == 1,:);
            % FILTER BY DIRECTION
            left=data_video_smoothfilt(data_video_smoothfilt(:,4)>90 & data_video_smoothfilt(:,4)<270,:);
            right=data_video_smoothfilt(data_video_smoothfilt(:,4)<90 | data_video_smoothfilt(:,4)>270,:);
            
            [SmoothRateMap_Right(i,:),num_spikes_Right]=RightVsLeft(right,spks_VEL(:,1),track_length,sampleRate);
            [SmoothRateMap_Left(i,:),num_spikes_Left]=RightVsLeft(left,spks_VEL(:,1),track_length,sampleRate);
            
            DirectionalityIndex=(num_spikes_Left-num_spikes_Right)/(num_spikes_Left+num_spikes_Right);
            Overall_DirectionalityIndex(i,:)=[num_spikes_Left num_spikes_Right];
        end  
        
%             % SPLIT UP LINEAR TRACK DATA BY DIRECTION
%     
%     if strcmp(linear_track,'yes')==1; 
%         mazeType=2; 
%         data_video_smoothfiltmaster=data_video_smoothfilt;
%     else
%         mazeType=1;
%     end
%     for maze=1:mazeType
%         if strcmp(linear_track,'yes')==1;
%             if maze==1
%                 direction='right';
%                 data_video_smoothfilt=data_video_smoothfiltmaster(data_video_smoothfiltmaster(:,4)>90 & data_video_smoothfiltmaster(:,4)<270,:);
%             elseif maze==2
%                 direction='left';
%                 data_video_smoothfilt=data_video_smoothfiltmaster(data_video_smoothfiltmaster(:,4)<90 | data_video_smoothfiltmaster(:,4)>270,:);
%             end
%         end
        
        % extract spike locations from smoothed and velocity filtered data
        spks_VEL = data_video_smoothfilt(data_video_smoothfilt(:,5) == 1,:);
        
        % filter by less than 5 spikes & delete old files
        if size(spks_VEL,1)<5
            clear spks_VEL % REMOVE SPIKES
            data_video_smoothfilt(:,5)=[]; % REMOVE SPIKES FROM VIDEO DATA
            counter = counter + 1;
            [filepath, filename] = fileparts(tfile{i});
            disp(['CELL ',filename, ' LESS THAN 5 SPIKES >>>>> DELETING ASSOCIATED FILES'])
            %             delete([filepath filesep filename '_spikeData.mat']);
            %             delete([filepath filesep filename '_spikeFigure.jpg']);
            %             delete([filepath filesep filename '_spikeAngle.jpg']);
            %             delete([filepath filesep filename '_PolarPlot.jpg']);
            %             % delete newer files
            %             delete([filepath filesep filename sprintf('S%d',event) '_spikeData.mat']);
            %             delete([filepath filesep filename sprintf('S%d',event) '_spikeFigure.jpg']);
            %             delete([filepath filesep filename sprintf('S%d',event) '_spikeAngle.jpg']);
            %             delete([filepath filesep filename sprintf('S%d',event) '_PolarPlot.jpg']);
            continue
        end
        %%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%

        % LFP ANALYSIS
        % Create name of CSC file from cell and path name
        [~, filename] = fileparts(tfile{i}); tetrode=regexp(filename,'.','match'); eegfile=cellstr(strcat(path,filesep,'CSC',tetrode(4),'.ncs'));
        
        % Compute Phase Lock
        [Stats,ThPrecess,ThetaStats]=EEGWorking2(eegfile{1},spks_VEL,StartofRec,EndofRec,event,linear_track,track_length);
        
        % BIN DATA
        [FilledRateMatrix,Occ1ZeroLogical,nBinsx,nBinsy,occ] = bindata( [data_video_smoothfilt(:,3),data_video_smoothfilt(:,2)],sampleRate,spks_VEL,linear_track,track_length);
        
        % smooth raw binned rate matrix using a guassian
        if isequal(linear_track, 'no')
            filtWidth = [3 3];
            filtSigma = 1;
            imageFilter=fspecial('gaussian',filtWidth,filtSigma);
            SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout');
        else
            filtWidth = [1,3];
            filtSigma = 1;
            imageFilter=fspecial('gaussian',filtWidth,filtSigma);
            SmoothRateMap = nanconv(FilledRateMatrix,imageFilter, 'nanout','1d');
            % SAVE ALL RATEMAP FOR POPULATION CODE
            SmoothRateMapAll(i,:)=SmoothRateMap;
        end
        % calculate information content
        rY = reshape(SmoothRateMap,nBinsx*nBinsy,1); % reshape data into column
        rY(isnan(rY))=0;
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
        DistanceFromTrackEnd = distanceMin*5;
        
        % COHERENCE
        if length(rY)>1; Coherence=corr2(FilledRateMatrix,SmoothRateMap); else Coherence=NaN; end
        
        % OVERALL FIRING RATE
        OverallFR = (length(spks_VEL)/length(data_video_smoothfilt))*sampleRate;
        
        % NUMBER OF SPIKES
        nSpikes=length(spks_VEL);
        
        % calculate sparsity
        [sparsity] = Sparsity(rY'); % from NSMA toolbox
        
        % calculate percent of active bins (estimate of field size - Royer et al., 2010)
        r_50perc = r_max*0.20;
        r_50percAll = find(rY > r_50perc);
        NumbActiveBins = numel(r_50percAll);
        FieldWidth = NumbActiveBins;
        
        % SPIKE DIRECTION (from Shawn Winter's code)
        [mean_vector_length,peak_Firing_Rate,preferred_Direction,halfPeak,Directional_Range_HalfWidth,Direct_infoContent,BinsNbSpikes,BinsAngle3,BinsAngle]...
            = HDCell(data_video_smoothfilt(:,5),data_video_smoothfilt(:,4),sampleRate);
        
        % -------------------------------CREATING PLOTS----------------------------
        if figures==1
            if isequal(linear_track, 'yes')
                % plot spike on path
                figure (i), subplot(5,1,1), plot(data_video_smoothfilt(:,2), data_video_smoothfilt(:,3), 'LineWidth', 1, 'color', 'k');
                hold on; scatter(spks_VEL(:,2), spks_VEL(:,3), 35, 'filled', 'r'); box off; axis off
                title(['Spike on Path /',' nSpikes: ',num2str(nSpikes)]);
                
                % plot filled and smoothed rate map
                figure (i), subplot(5,1,2), h = pcolor([SmoothRateMap;SmoothRateMap]);
                hold on; shading interp; colormap jet; axis off
                hold on; box off; set(h, 'EdgeColor', 'none'); axis image
                title(['Smoothed RateMap',' InfoContent: ',num2str(InformationContent)]);
                
                
                figure (i), subplot(5,1,3), h = area(SmoothRateMap(1,:),'LineWidth',2,'EdgeColor',[0,0,0]+0.4,'FaceColor', [0,0,0]+0.8);
                box off; xlim([1 nBinsx]);
                title(['Rate Plot',', Peak Rate: ',num2str(PeakRate)]);
                set(figure (i),'Position',[842 42 838 924]);
                
                % PLOT SCATTERED PHASE PRECESSION
                trackX=ThPrecess.scatteredPH(:,1).*pixelDist;
                figure(i); subplot(5,1,4); plot(trackX, ThPrecess.scatteredPH(:,2), 'k.')
                hold on; plot(trackX, ThPrecess.scatteredPH(:,2)+360, 'r.')
                ylim([0 720]); xlim([min(trackX),max(trackX)])
                set(gca, 'YTick', [0;240;480;720],'Box','off');
                title(['PhasePrecession',' R-Squared: ',num2str(ThPrecess.RSquared),...
                    ', Slope: ',num2str(ThPrecess.slope),', Corr: ',num2str(ThPrecess.Correlation)])
                hold off
                
                % PLOT SMOOTHED RATEMAP
                figure(i); subplot(5,1,5); h=pcolor(ThPrecess.smoothedPHmap);
                shading interp; colormap jet; hold on; box off; axis off; set(h, 'EdgeColor', 'none');
                title(['Mean Firing Rate: ',num2str(ThPrecess.meanFR),', DOM: ',num2str(ThPrecess.DOM)])
                
                % subplot5 = subplot(5,1,5,'Parent',figure(i),'YTickLabel',{'0','240','480','720'},'YTick',[0 240 480 720]);
%                 surface('Parent',subplot5,'AlignVertexCenters','on','FaceColor','interp','EdgeColor','none',...
%                     'CData',expand(ThPrecess.smoothedPHmap,[720/size(ThPrecess.smoothedPHmap,1),1]),...
%                     'ZData',expand(ThPrecess.smoothedPHmap,[720/size(ThPrecess.smoothedPHmap,1),1]),...
%                     'YData',1:720,...
%                     'XData',1:size(ThPrecess.smoothedPHmap,2));
%                 colormap jet; ylim([0 720]); xlim([1 size(ThPrecess.smoothedPHmap,2)]); shading interp;
%                 title(['Mean Firing Rate: ',num2str(ThPrecess.meanFR),', DOM: ',num2str(ThPrecess.DOM)])
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
                shading interp; colormap jet; axis off; hold on; box off; set(h, 'EdgeColor', 'none'); axis image
                title(['Smoothed Rate Map, IC: ',num2str(InformationContent)]);
                
                %smooth spike data
                [pdf, sigma] = circ_ksdensity(spks_VEL(:,4), 0:359,'msni'); 
                
                % Plot Smoothed firing rate x HEAD DIRECTION
%                 figure (i+2), plot(BinsAngle3,BinsNbSpikes,'LineWidth',2,'color','k') %plot(BinsAngle3,BinsNbSpikes,'LineWidth',2,'color','k')
                figure (i),subplot(2,2,3); plot(pdf,'LineWidth',2,'color','k') 
                axis tight; hold on; xlim ([0 360]); box off
                title(['Tuning Curve, Peak Rate: ',num2str(peak_Firing_Rate),' D_IC: ',num2str(Direct_infoContent)])
                
                % Firing rate x HD polar plot for the nonsmoothed data above
                figure(i); subplot(2,2,4); polarplot = polar(BinsAngle([1:60 1]),BinsNbSpikes([1:60 1]),'b');
                set(polarplot, 'linewidth',3,'color','k'); axis off
                title(['Polor Plot, MeanVecLength: ',num2str(mean_vector_length),' Pref_Dir: ',num2str(preferred_Direction)]);
                set(0,'Showhiddenhandles','on')
                
                % ---------CODE FOR PUBLICATION FIGURE--------
                % extrastuff = setdiff(get(gca,'children'),polarplot);
                % delete(extrastuff)
                % horizontal=line([-max(BinsNbSpikes) max(BinsNbSpikes)],[0 0]);
                % vertical=line([0 0],[-max(BinsNbSpikes) max(BinsNbSpikes)]);
                % set(horizontal,'linewidth',2,'color','k');
                % set(vertical,'linewidth',2,'color','k');
                %---------------------------------------------
                set(figure(i),'Position',[686 325 977 619]);
                fig1 = figure(i);
            end
        end
        
        % SAVING JPGS AND .MAT FILES TO TT FOLDER
        cd(mclustpath);
        warning('off', 'MATLAB:Figure:FigureSavedToMATFile');
        
        if exist('timestamps','var'); % SAVE WITH S# IF EVENTS EXIST
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
        
        keep('path','data_video_smoothfilt','SpikeFile','mclustpath','tfile','track_length',...
            'linear_track','counter','SmoothRateMap_Right','SmoothRateMap_Left'...
            ,'figures','S','sampleRate','i','Overall_DirectionalityIndex','StartofRec',...
            'EndofRec','event','timestamps','vel_cmPerSec','table2','data_video_smoothfilt2',...
            'lowercut','highercut','SmoothRateMapAll','pixelDist');
        
        data_video_smoothfilt(:,5)=[]; % REMOVE SPIKES FROM VIDEO DATA
        counter = counter + 1;
        close all
    end
%     end % put at end of script for linear track direction

    % DIRECTIONAL MATRIX CREATION
    if isequal(linear_track, 'yes') && exist('SmoothRateMap_Right','var')==1
        if exist('timestamps','var')==1; timestampz=1; else timestampz=0; end
        Overall_DirectionalityIndex=DirectionalMatrix(tfile{1},event,timestampz,SmoothRateMap_Right,SmoothRateMap_Left,Overall_DirectionalityIndex);
        DirectionalMatrix(tfile{1},event,timestampz,SmoothRateMapAll);
    end
end
disp 'DONE'




