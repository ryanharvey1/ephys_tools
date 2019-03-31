function  postprocessMClust(  path, track_length, linear_track )
%%%%%%%% DATA EXTRACTION & PREPROCESSING from Neuralynx Files %%%%%%%%%%%%%
% Can handle circular arena and linear track data
% Ryan H and Ben C May 2016

path = 'D:\2016-08-06_12-34-42';
track_length = 120;
linear_track = 'yes';

% load path / video file / spike files

% add needed functions to path 
try
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\spikeCode'));  
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\Cell analysis'));  
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\chronux_2_11'));  
catch
end
VTfile = FindFiles('*.nvt', 'StartingDirectory', path); % for Nlx2MatVT compiler

% ReadVideoTrackingFile3(path); % extracts raw xy for each led, timestamps(might be incorrect) and saves txt output to parent folder
% VTfile = FindFiles('VT1.txt', 'StartingDirectory',path);  % locates .txt
% table=table2array(readtable('VT1.txt')); % reads table and converts it into an array
% 
% % Not best way to do this, but follow my logic below
% readnvt=ReadNVT('VT1.nvt'); % extracts xy locations,timestamps, angles
% ExtractedAngle=readnvt.Angle; % specifies angle from ReadNVT function
% VTTimeStamps = readnvt.TimeStamp; % specifies timestamp from ReadNVT function because timestamp from ReadVideoTrackingFile3 error 
% table(:,1) = []; % removes old timestamp column
% ExtractedX = readnvt.Xloc;
% ExtractedY = readnvt.Yloc;
% VTfile=[VTTimeStamps,table,ExtractedX,ExtractedY,ExtractedAngle]; %combines the timestamp, raw xy for each led with the respective angle of rat's head
% 
% c = find(VTfile(:,2)>0); c = VTfile(c,:); % removes 0 from led xy
% c = find(VTfile(:,3)>0); c = VTfile(c,:); % faster if used logicals instead of find (update in the future)
% c = find(VTfile(:,4)>0); c = VTfile(c,:);
% c = find(VTfile(:,5)>0); c = VTfile(c,:);
% % We not have variable c which contains timestamp, xy for each led, and angle of rats head 
% c(:,[2,3,4,5]) = []; % removes raw xy for both leds
% c=[VTTimeStamps,ExtractedX,ExtractedY,ExtractedAngle];

try
    SpikeFile=FindFiles('*SS*.txt','StartingDirectory',path); % For Spikesort 3D
catch
end
if exist('SpikeFile')==2; % For Spikesort 3D
    SS3D=2;
end
if  exist('SpikeFile')==1;
%     try 
    addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox')); % add NSMA to path
%     catch  
%     try 
%     NSMA_Toolbox = ('/Users/RyanHarvey/Dropbox/MATLAB/NSMA_Toolbox'); % add NSMA to path
%     addpath(genpath(NSMA_Toolbox));
%     catch 
%     end
%     try
    mclustpath=cd([path,'p',filesep 'TT']); % sets path to p file & looks within TT folder

    ls '*.t' % disp .t files
    tfile=FindFiles('*.t');

for i=1:length(tfile); % Counts # of .t files (# of clusters) and sets iteration length
    SS3D=1; % For MClust
    addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox')); % add NSMA to path
    S =LoadSpikes(tfile); % loads .t into ts cells
    t = Data(S{i}); % converts ts cells to numerals 
    SpikeFile=100*t; % conversion from MClust output to our current code 
    rmpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox')); % Remove NSMA from path (different Nlynx code)

% for i = 1:length(SpikeFile);  % For Spikesort 3D (commented out, still working on it.)

% video frame rate in Hz
sampleRate = 30;

% extract video data with mex
% Nlx2MatVT_v3 for mac compiler *** Nlx2MatVT for PC compiler 
try
[VTTimeStamps, ExtractedX, ExtractedY, ExtractedAngle] = Nlx2MatVT(VTfile{1}, [1 1 1 1 0 0], 0, 1);  % for PC
catch 
end
try
[VTTimeStamps, ExtractedX, ExtractedY, ExtractedAngle] = Nlx2MatVT_v3(VTfile{1}, [1 1 1 1 0 0], 0, 1);  % for mac                                                    
catch
end
% [VTTimeStamps, ExtractedX, ExtractedY, ExtractedAngle]=c; % bypass nlx2matvt compiler

ts_video = VTTimeStamps'; 
angle_video = ExtractedAngle';
data_video = [ts_video, ExtractedX', ExtractedY', angle_video]; % Graph this angle data / smooth it

%Filter out tracking issues
data_video1 = (data_video(:,2)>1); % can be changed based on maze config
data_video = data_video(data_video1,:);
data_video1 = (data_video(:,3)>1); % can be changed based on maze config
data_video = data_video(data_video1,:);

% Smoothing data
x_smooth = runline(data_video(:,2),10,1); 
y_smooth = runline(data_video(:,3),10,1);
data_video_smoothfilt = [data_video(:,1) x_smooth y_smooth data_video(:,4)];

%Interpolate 
% if SS3D==2
%     s1 = importdata(SpikeFile{i});       % for spikesort 3D text files
% end
if SS3D==1
    s1 = SpikeFile;                    % for MClust .t files
end
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


%%%%%%%%%%%% DATA ANALYSIS %%%%%%%%%%%%%%%%%%%%%

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
iVel = find(vel_cmPerSec >= 0); % 5cm/sec is used by Stensola et al 2012
data_video_smoothvelfilt = data_video_smoothfilt(iVel,:);

% extract spike locations from smoothed and velocity filtered data
spksi = find(data_video_smoothvelfilt(:,5) == 1); 
spks_VEL = data_video_smoothvelfilt(spksi,:);

% create array for binning data
if isequal(linear_track, 'yes')
    nBinsx = round(track_length/4); %~4x4cm bin size
    % nBinsx = 30; % ~4x4cm assuming a ~120cm track
    nBinsy = 2;
elseif isequal(linear_track,'no')
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

% calculate sparsity
[sparsity] = Sparsity(rY'); % from NSMA toolbox

% calculate percent of active bins (estimate of field size - Royer et al., 2010)
r_50perc = r_max*0.50;
r_50percAll = find(rY > r_50perc);
NumbActiveBins = numel(r_50percAll);
FieldWidth = NumbActiveBins;

% Remove unsampled bins for plotting
SmoothRateMap(Occ1ZeroLogical) = NaN;
FilledRateMatrix(Occ1ZeroLogical) = NaN;

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
% % images / data to retain
% if SS3D==2
% %     cd(path);
%     [filepath, filename] = fileparts(SpikeFile{i});
%     save([filepath filesep filename '_spikeData.mat']);
%     print(figure (i), '-djpeg', filename);
%     close(figure (i)); 
% end

if SS3D==1
    cd(mclustpath);
    [filepath, filename] = fileparts(tfile{i});
    saveas(fig1,[filepath filesep filename '_spikeFigure.jpg']);
    if isequal(linear_track, 'no')
    saveas(fig2,[filepath filesep filename '_RateMapFigure.jpg']);
    end
    save([filepath filesep filename '_spikeData.mat']);
%     close all
%     cd([mclustpath,'p','\TT']);
%     print(figure (i), '-djpeg', filename);
%     close(figure (i)); 
    keep('path','VTfile','SpikeFile','mclustpath','tfile','track_length','linear_track'); 
close all
end
end
end
% end
% end





