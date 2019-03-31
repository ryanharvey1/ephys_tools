%% load video file
path = '/Users/ryanharvey/Dropbox/TestData';
VTfile = FindFiles('*.nvt', 'StartingDirectory', path);
SpikeFile = FindFiles('TT4 RH_SS_01.txt', 'StartingDirectory', path);

% video frame rate in Hz
sampleRate = 30;

% extract video data with mex
% switched Nlx2MatVT to Nlx2MatVT_v3 for mac compiler 
[VTTimeStamps, ExtractedX, ExtractedY, ExtractedAngle] = Nlx2MatVT_v3(VTfile{1}, [1 1 1 1 0 0], 0, 1);                                                      
xy_pos = [ExtractedX', ExtractedY'];
ts_video = VTTimeStamps';
angle_video = ExtractedAngle';
data_video = [ts_video, ExtractedX', ExtractedY', angle_video];

%Filter out tracking issues
data_video1=(data_video(:,2)>10); 
data_video=data_video(data_video1,:);
data_video1=(data_video(:,3)>100); 
data_video=data_video(data_video1,:);

% Smoothing data
x_smooth=runline(data_video(:,2),7,1);
y_smooth=runline(data_video(:,3),7,1);
data_video_smoothfilt = [data_video(:,1) x_smooth y_smooth data_video(:,4)];
% 
% % velocity of rat from smoothed xy data
% vel_x = diff(x_smooth); % vel units are pixels/frame
% vel_y = diff(y_smooth);
% vel_abs = sqrt(vel_x.^2 + vel_y.^2); % scalar length of velocity vector = "scalar velocity" in pixels/frame
% vel_cmPerSec = vel_abs * 2.1 * sampleRate; % ~2.1 cm/pixel for Taube lab according to Shawn W  
% 
% a = vel_cmPerSec(1,:);
% vel_cmPerSec = [a; vel_cmPerSec];
% data_video_smoothfilt = [data_video(:,1) x_smooth y_smooth data_video(:,4) vel_cmPerSec];
% 
% % find xy coords and spikes when rat is above velocity threshold
% iVel = find(data_video_smoothfilt(:,5) >= 0); % 5cm/sec is used by Stensola et al 2012
% % spksf = data_video_smoothfilt(:,4); % **************Prob not right!********
% % dataVEL = [x_smooth y_smooth spksf];
% data_video_smoothvelfilt = data_video_smoothfilt(iVel,:);

% % extract spike locations from smoothed and velocity filtered data
% spksi = find(dataVEL(:,3) == 1); 
% spks_VEL = dataVEL(spksi,:);

%Interpolate 
s1 = importdata(SpikeFile{1});
ts_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,1), s1, 'nearest');
x_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,2), s1, 'nearest');
y_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,3), s1, 'nearest');
angle_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,4), s1, 'nearest');
 
% create array for binning data
nBins = 15; % adjusted so the bins are ~5x5cm as in Koenig et al 2011
MinY = min(data_video_smoothfilt(:,3));
MaxY = max(data_video_smoothfilt(:,3));
MinX = min(data_video_smoothfilt(:,2));
MaxX = max(data_video_smoothfilt(:,2));
edges{1} = linspace(MinY, MaxY, nBins+1);
edges{2} = linspace(MinX, MaxX, nBins+1);

% bin occupancy data
occMatrix = [data_video_smoothfilt(:,3),data_video_smoothfilt(:,2)];
Omatrix = hist3(occMatrix,'Edges',edges);
Omatrix(1,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/sampleRate;

% NEW find pixels with zero occupancy- these will get blanked out out in the
% rate maps
Occ1NonZeroLogical = Omatrix~=0;
Occ1ZeroLogical = Omatrix == 0;

% bin spike data
spikeMatrix = [y_spike, x_spike];
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

% NEW
SmoothRateMap(Occ1ZeroLogical) = NaN;
FilledRateMatrix(Occ1ZeroLogical) = NaN;

% plot filtered and smoothed rate map 
figure (1), subplot(2,1,1), plot(data_video_smoothfilt(:,2), data_video_smoothfilt(:,3), 'LineWidth', 2, 'color', [0,0,0]+0.8);
hold on;
scatter(x_spike, y_spike, 45, 'filled', 'k');
box off
axis off
title('Spike (black dots) on Path (gray)');

% plot filled and smoothed rate map 
NaNsForSmoothedR = NaN(1,nBins+1);
NaNsForSmoothedC = NaN(nBins,1);
SmoothRateMapGraph = [SmoothRateMap,NaNsForSmoothedC;NaNsForSmoothedR];
figure (1), subplot(2,1,2), h = pcolor(SmoothRateMapGraph);
colormap jet
axis off
hold on
% colorbar
box off
set(h, 'EdgeColor', 'none');
% set(gca,'YDir','reverse');
title('Smoothed Rate Map');

% [filepath, filename] = fileparts(SpikeFile);
% savefig([filepath filesep filename '.fig'])
% save([filepath filesep filename '.mat']);

%% 
clear all
close all
