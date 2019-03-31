
%%%%%%%% DATA EXTRACTION & POSTPROCESSING from Neuralynx Files %%%%%%%%%%%%%
% postprocessMClust_v5_batch
% *THIS VERSION CYCLES THROUGH ALL DATA*
% THIS SCRIPT CREATES FIGURES AND STATS FOR PLACE & HD CELL DATA
% CAN HANDLE CIRCULAR ARENA & LINEAR TRACK DATA 
%
% OUTPUTS: RATE MAP / SPIKE ON PATH / TUNING CURVE / POLAR PLOT
% Ryan Harvey & Ben Clark 2016  
clear, clc , close all
%#######################################################
linear_track = 'yes'; % <<<<<<<<<<<<<<< 'yes' or 'no'
%#######################################################

rats={'RH11','RH13','RH14','RH16'}; % LIST RATS YOU WANT TO RUN

for irats=1:length(rats)
    if ismac==1
        parent = strcat('/Users/RyanHarvey/OneDrive - University of New Mexico/Test_Sessions/',rats(irats)); % CHANGE BASED ON PARENT FOLDER
    else
        parent = strcat('D:\Place_Cell_Data\PAE_Rat\',rats(irats)); % CHANGE BASED ON PARENT FOLDER
    end
disp(['CYCLING THROUGH RAT:',char(rats(irats))])
parent=char(parent);
structdir=dir(parent);
for I=1:length(structdir) % 1 TO # OF FILES IN DIR
    if structdir(I).isdir && structdir(I).name(1) ~= '.' && any(regexp(structdir(I).name,'p$'))~=1; % IF NOT '.' & NOT P FILE
        cd([parent filesep structdir(I).name,'p',filesep 'TT']); % CD TO TT FOLDER       
            if any(size(dir([parent filesep structdir(I).name,'p',filesep 'TT' filesep '*.t' ]),1))==1 % CHECK FOR MCLUST FILES
                path=[parent filesep structdir(I).name]; % SET PATH TO DATA THAT HAS BEEN THROUGH MCLUST
            else
                continue
            end
    else
        continue
    end
    
%% SPECIFY TRACK LENGTH >>> 90 (pre 7/14/16) or 120 (7/14/16 and on)
SplitPath=strsplit(path,'-'); % SPLIT PATH UP BY -
SplitPath2=strsplit(char(SplitPath(3)),'_'); % SPLIT PATH UP BY _
if str2double(SplitPath(2))<=7; % IF MONTH IS 7 OR LESS
    track_length = 90;
    if str2double(SplitPath(2))==7;
        if str2double(SplitPath2(1))<13; % IF DAY IS LESS THAN 14
            track_length = 90;
        else str2double(SplitPath2(1))>13;
            track_length = 120; 
        end
    end
else str2double(SplitPath(2))>7;
    track_length = 120; 
end

%% ADD TOOLBOXES TO PATH
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\spikeCode'));  % PC
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\Cell analysis'));  
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\chronux_2_11'));  
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\PAE_Project')); 
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\CircStat2012a'));


addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/spikeCode'));  % MAC
addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/Cell analysis'));  
addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/chronux_2_11'));  
addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/PAE_Project')); 
addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/CircStat2012a'));
  
%% BYPASSING Nlx2MatVT COMPILER WITH ReadVideoTrackingFile3 (TAKES A WHILE TO RUN)
ReadVideoTrackingFile3(path); % extracts raw xy for each led
VTfile = FindFiles('VT1.txt', 'StartingDirectory',path);  % locates newly created VT1.txt file
table=textread(cell2mat(VTfile)); % reads table 

% BYPASSING Nlx2MatVT COMPILER WITH READNVT
readnvt=ReadNVT('VT1.nvt'); % extracts xy locations,timestamps, angles
VTTimeStamps = readnvt.TimeStamp; % specifies timestamp from ReadNVT function because timestamp from ReadVideoTrackingFile3 error 
table(:,1) = []; % removes old timestamp column
ExtractedX = readnvt.Xloc;
ExtractedY = readnvt.Yloc;
ExtractedAngle = atan2d(table(:,4)-table(:,2),table(:,3)-table(:,1)) + 360*((table(:,4)-table(:,2))<0); % GETS ANGLE FROM RAW XY FOR EACH LED
VTfile=double([VTTimeStamps,table,ExtractedX,ExtractedY,ExtractedAngle]); %combines the timestamp, raw xy for each led with the respective angle of rat's head

%% IF TRACKING NON-DETECT>>>REPLACE ANGLE WITH NAN
for izeroes=1:length(VTfile)
    if VTfile(izeroes,2)<=0
        VTfile(izeroes,8)=NaN; 
    end
    if VTfile(izeroes,4)<=0
        VTfile(izeroes,8)=NaN;
    end
end

% REMOVES TABLE FROM VTFILE
VTfile(:,[2,3,4,5])=[]; 

% CALC NUMBER AND PERCENT OF NON-DETECTS
SumofNaN=sum(isnan(VTfile(:,4)));
PercentofNaN=(SumofNaN/length(VTfile))*100;

% IF *FIRST* ANGLE IS NAN, USE AVERAGE OF FIRST ANGLES
if isnan(VTfile(1,4))==1
    for inan=1:length(VTfile)
        if isnan(VTfile(inan+1,4))==0 % COUNTS # OF ROWS OF NANs
            n=inan+1; % # OF NANS
            break
        end
    end
    VTfile(1,4)=nanmean(VTfile(1:n,4)); % AVG OF n
end
% IF *LAST* ANGLE IS NAN, USE AVERAGE OF LAST ANGLES
if isnan(VTfile(end,4))==1
    for inan=length(VTfile):-1:1
        if isnan(VTfile(inan-1,4))==0 % COUNTS # OF ROWS OF NANs
            n=inan-1; % # OF NANS
            break
        end
    end
    VTfile(end,4)=nanmean(VTfile(n:end,4)); % AVG OF n
end

% FILLING IN NaNs
[VTfile]=FillNaN2(VTfile,4); 
disp(['Filling in ', num2str(SumofNaN), ' non-detects, which is ', num2str(PercentofNaN), ' percent of your data.'])

% IF NaNs ARE STILL FOUND (they shouldn't be), INDICATE THE ROWS WHERE THEY ARE
if isnan(VTfile)==1
    for inan=1:length(VTfile)
        if isnan(VTfile(inan,4))==1
            fprintf('NaN found on row %d \n',inan);
        end
    end
end
%% LOCATE T FILES
mclustpath=cd([path,'p',filesep 'TT']); % sets path to p file & looks within TT folder for tfiles
    
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox')); % add NSMA to path for FindFiles
addpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox'));
    
ls '*.t' % disp .t files
tfile=FindFiles('*.t');
   
%% CYCLE THROUGH CLUSTERS AND OUTPUT FIGURES + .MAT FILES
counter = 1;
for i=1:length(tfile); % Counts # of .t files (# of clusters) and sets iteration length
    addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox')); % add NSMA to path
    addpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox')); % add NSMA to path
    S =LoadSpikes(tfile); % loads .t into ts cells
    t = Data(S{i}); % converts ts cells to numerals 
    SpikeFile=100*t; % conversion from MClust output to our current code 
    
    rmpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)
    rmpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)

% video frame rate in Hz
sampleRate = 30; % 30 TIMESTAMPS = 1 SECOND

%Filter out tracking issues
VTfile = VTfile(VTfile(:,2)>1,:); % ONLY KEEP TRACKING DATA ABOVE 1 
VTfile = VTfile(VTfile(:,3)>1,:); 

% Smoothing data
x_smooth = runline(VTfile(:,2),10,1); 
y_smooth = runline(VTfile(:,3),10,1);
data_video_smoothfilt = [VTfile(:,1) x_smooth y_smooth VTfile(:,4)];

%Interpolate 
s1 = SpikeFile; % for MClust .t files

ts_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,1), s1, 'nearest');
% x_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,2), s1, 'nearest');
% y_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,3), s1, 'nearest');
angle_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,4), s1, 'nearest'); 

%Create a vector of 0's and 1's (1 = ts in which spike occured) and add to
%array of data
[~,ia,~] = intersect(data_video_smoothfilt(:,1), ts_spike);
N = size(data_video_smoothfilt(:,1));
spikeVec = zeros(N(1,1),1);
spikeVec(ia) = 1;
data_video_smoothfilt = [data_video_smoothfilt spikeVec(:,1)];


%% %%%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%

% velocity of rat from smoothed xy data
vel_x = diff(data_video_smoothfilt(:,2)); % vel units are pixels/frame
vel_y = diff(data_video_smoothfilt(:,3));
vel_abs = sqrt(vel_x.^2 + vel_y.^2); % scalar length of velocity vector = "scalar velocity" in pixels/frame
xmin = min(vel_x);
xmax = max(vel_x);
xdist = xmax - xmin;
if isequal(linear_track,'no')
    track_length = 24;
end
pixelDist = track_length/xdist;
vel_cmPerSec = vel_abs * pixelDist * sampleRate; % ~2.1 cm/pixel for Taube lab according to Shawn W  

% find xy coords and spikes when rat is above velocity threshold
data_video_smoothvelfilt = data_video_smoothfilt((vel_cmPerSec >= 0),:); % 5cm/sec is used by Stensola et al 2012

% extract spike locations from smoothed and velocity filtered data
spks_VEL = data_video_smoothvelfilt(data_video_smoothvelfilt(:,5) == 1,:);

% create array for binning data
if isequal(linear_track, 'yes')
    nBinsx = round(track_length/4); %~4x4cm bin size
    % nBinsx = 30; % ~4x4cm assuming a ~120cm track
    nBinsy = 2;
elseif isequal(linear_track,'no') % FOR ARENA
    nBinsx = 24;
    nBinsy = 24;
end

MinY = min(data_video_smoothvelfilt(:,3));
MaxY = max(data_video_smoothvelfilt(:,3));
MinX = min(data_video_smoothvelfilt(:,2));
MaxX = max(data_video_smoothvelfilt(:,2));
edges{1} = linspace(MinY, MaxY, nBinsy+1);
edges{2} = linspace(MinX, MaxX, nBinsx+1);

% bin occupancy data
occMatrix = [data_video_smoothvelfilt(:,3),data_video_smoothvelfilt(:,2)];
Omatrix = hist3(occMatrix,'Edges',edges);
Omatrix(1,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/sampleRate;

% find pixels with zero occupancy- these will get blanked out out in the
% rate maps
Occ1NonZeroLogical = Omatrix~=0;
Occ1ZeroLogical = Omatrix == 0;

% bin spike data
spikeMatrix = [spks_VEL(:,3), spks_VEL(:,2)];
Smatrix = hist3(spikeMatrix,'Edges',edges);
Smatrix(1,:) = [];
Smatrix(:,end) = [];

% divide binned spikes by occupancy to get rate maps and store data
% (can use occ+eps instead of removing 0's)
BFiringRateMatrix = Smatrix./occ;
FilledRateMatrix = BFiringRateMatrix;
FilledRateMatrix(isnan(FilledRateMatrix)) = 0;
FilledRateMatrix(isinf(FilledRateMatrix)) = 0;

% smooth raw rate histogram using a guassian (from NSMA toolbox)
[SmoothRateMap] = smooth(FilledRateMatrix,1,5,1,5); % smooth with 5x5pt window and 1pt step

% calculate information content
rY = reshape(SmoothRateMap,nBinsx*nBinsy,1); % reshape data into column
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
DistanceFromTrackEnd = distanceMin*4;

% OVERALL FIRING RATE
OverallFR = (length(s1)/length(VTfile))*sampleRate;

% calculate sparsity
[sparsity] = Sparsity(rY'); % from NSMA toolbox

% calculate percent of active bins (estimate of field size - Royer et al., 2010)
% r_50perc = r_max*0.50;
r_50perc = r_max*0.25; % CHANGED TO 25 PERCENT (ACTIVE FIELD IS BIN WITH 25% OF THE MAXIMUM FIRING RATE) RYAN H 10/27/16
r_50percAll = find(rY > r_50perc);
NumbActiveBins = numel(r_50percAll);
FieldWidth = NumbActiveBins;

% Remove unsampled bins for plotting
SmoothRateMap(Occ1ZeroLogical) = NaN;
FilledRateMatrix(Occ1ZeroLogical) = NaN;

%% SPIKE DIRECTION 

% extract spike directions from smoothed and velocity filtered data
spks = data_video_smoothvelfilt(:,5);
angle = data_video_smoothvelfilt(:,4);

% number of directional bins
dBins = 60;
        
% Binning of directional data can also be done using histc function
    for j = 1:dBins;
        k = 0:dBins;
        % number of instances where the head of the rat was found in each of the 60 possible orientations
        ListOrientation = find(angle >= k(j)*6 & angle < (k(j)+1)*6);
            if length(ListOrientation) < 1; % if the number is less than 1 (i.e., 0), it is hard set to 1 to avoid division by 0
                nOrientation(j) = 1; 
                nSpikesOrientation(j) = 0; % 0 is assigned to the number of spikes for this orientation
            else
                nSpikesOrientation(j) = sum(spks(ListOrientation));
                nOrientation(j) = length(ListOrientation);
            end
    end
    
% transformed values from 1/60th of a sec to seconds
nOrientation2 = nOrientation./sampleRate; % 30ms for the time in s

% calculates the spikes/sec for each 6° bin
BinsNbSpikes = nSpikesOrientation./nOrientation2;

% trace the graph of the discharge rate / direction to the cells 1
BinsAngle = 0.052358333:0.104716667:6.283;
BinsAngle3 = 3:6:360;

% calculate head direction cell properties
mean_vector_length = circ_r(BinsAngle',BinsNbSpikes',circ_ang2rad(6)); % mean vector length based on binned firing rates
peak_Firing_Rate = max(BinsNbSpikes); % peak firing rate
pfdi = find(BinsNbSpikes(1,:) == peak_Firing_Rate); 
preferred_Direction = BinsAngle3(pfdi); % preferred firing direction
halfPeak = peak_Firing_Rate/2;
hpi = find(BinsNbSpikes(1,:) >= halfPeak); 
Directional_Range_HalfWidth_bins = BinsAngle3(hpi);
Directional_Range_HalfWidth = max(Directional_Range_HalfWidth_bins) - min(Directional_Range_HalfWidth_bins);

% calculate directional information content (from Taube & Muller 1998)
probOrient = nOrientation2./sum(nOrientation2); % probability of occupancy
overall_rate = (sum(nSpikesOrientation))/(sum(nOrientation2));
reIC = BinsNbSpikes'./overall_rate;
log_IC = log2(reIC);
ij= find(isinf(log_IC)); % find -Inf's (log(0)) and replace with 0's (based on code from McN lab)
log_IC(ij) = 0;
ICi = probOrient'.*reIC.*log_IC;
Direct_infoContent = sum(ICi);

%% -------------------------------CREATING PLOTS----------------------------
if isequal(linear_track, 'yes')
% plot spike on path
figure (i), subplot(3,1,1), plot(data_video_smoothvelfilt(:,2), data_video_smoothvelfilt(:,3), 'LineWidth', 2, 'color', [0,0,0]+0.8);
hold on;
scatter(spks_VEL(:,2), spks_VEL(:,3), 45, 'filled', 'k');
box off
axis off
title('Spike (black dots) on Path (gray)');

% plot filled and smoothed rate map 
% NaNsForSmoothedR = NaN(1,nBins+1);
% NaNsForSmoothedC = NaN(nBins,1);
% SmoothRateMapGraph = [SmoothRateMap,NaNsForSmoothedC;NaNsForSmoothedR];
figure (i), subplot(3,1,2), h = pcolor(SmoothRateMap);
colormap jet
axis off
hold on
% colorbar
box off
set(h, 'EdgeColor', 'none');
% set(gca,'YDir','reverse');
title('Smoothed Rate Map');
axis image

figure (i), subplot(3,1,3), h = area(SmoothRateMap(1,:),'LineWidth',2,'EdgeColor',[0,0,0]+0.4,'FaceColor', [0,0,0]+0.8);
box off
xlim([1 nBinsx]);
title('Rate Plot');
fig1 = figure(i);
end 

%% FOR CIRCULAR ARENA
if isequal(linear_track, 'no')
% plot spike on path
figure (i), plot(data_video_smoothvelfilt(:,2), data_video_smoothvelfilt(:,3), 'LineWidth', 2, 'color', [0,0,0]+0.8);
hold on;
axis off
scatter(spks_VEL(:,2), spks_VEL(:,3), 45, 'filled', 'k');
box off
title('Spike (black dots) on Path (gray)');
axis image
fig1 = figure(i);

figure (i+1), h = pcolor(SmoothRateMap);
colormap jet
axis off
hold on
% colorbar
box off
set(h, 'EdgeColor', 'none');
% set(gca,'YDir','reverse');
title('Smoothed Rate Map');
axis image
fig2 = figure(i+1);
end

%% Plot firing rate x HEAD DIRECTION
figure (i+2), plot(BinsAngle3,BinsNbSpikes,'LineWidth',2,'color','k')
axis tight;
hold on
xlim ([0 360]);
box off
fig3=figure (i+2);

% Firing rate x HD polar plot for the data above
figure(i+3)
polarplot = polar(BinsAngle([1:60 1]),BinsNbSpikes([1:60 1]),'b');
set(polarplot, 'linewidth',3,'color','k');
axis off
set(0,'Showhiddenhandles','on')
% ---------CODE FOR PUBLICATION FIGURE--------
% extrastuff = setdiff(get(gca,'children'),polarplot);
% delete(extrastuff)
% horizontal=line([-max(BinsNbSpikes) max(BinsNbSpikes)],[0 0]);
% vertical=line([0 0],[-max(BinsNbSpikes) max(BinsNbSpikes)]);
% set(horizontal,'linewidth',2,'color','k');
% set(vertical,'linewidth',2,'color','k');
%---------------------------------------------
% title(['Peak Firing Rate:',num2str(peak_Firing_Rate),'Mean Vector Length',num2str(mean_vector_length),'Information Content',num2str(Direct_infoContent)])
fig4 = figure(i+3);

%% SAVING JPGS AND .MAT FILES TO TT FOLDER
    cd(mclustpath);
    warning('off', 'MATLAB:Figure:FigureSavedToMATFile');
    [filepath, filename] = fileparts(tfile{i});
    saveas(fig1,[filepath filesep filename '_spikeFigure.jpg']);
    if isequal(linear_track, 'no')
        saveas(fig2,[filepath filesep filename '_RateMapFigure.jpg']);
    end
    saveas(fig3,[filepath filesep filename '_spikeAngle.jpg']);
    saveas(fig4,[filepath filesep filename '_PolarPlot.jpg']);
    save([filepath filesep filename '_spikeData.mat']);
    keep('path','VTfile','SpikeFile','mclustpath','tfile','track_length','linear_track','counter','structdir','parent','rats','irats'); 
close all
fprintf('Just finished iteration #%d of %d \n', counter,length(tfile));
counter = counter + 1;
end
end
end
disp('DONE')






