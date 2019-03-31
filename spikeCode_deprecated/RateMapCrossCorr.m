% RateMapCrossCorr
%
% Input: 
%       - raw data output from Labview
%
% Output:
% 
% Created by Ben C September 2014

% identify path to data files
path = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__RateMapData/__control_data_26x26/mec';
ReadData(1,:) = FindFiles('bc217_s32_TT2_u2.txt', 'StartingDirectory', path);
ReadData(2,:) = FindFiles('bc217_s34_TT2_u2.txt', 'StartingDirectory', path);

% video frame rate in Hz
sampleRate = 60;

% maze parameters for non-detects
xmin = 0;
xmax = 255;
ymin = 0;
ymax = 255;

% min and max LED distance parameters
minLED = 4;
maxLED = 39;


% for i = 1:length(ReadData)
% load xy data
data = importdata(ReadData{1});
% extract coords, spikes, angle, and direction from ReadData output
rx = data.data(:,2); % red x-coords
ry = data.data(:,3); % red y-coords
gx = data.data(:,4); % green x-coords
gy = data.data(:,5); % green y-coords
spks = data.data(:,6); % spikes
angle = data.data(:,10); % angle in radians
distance = data.data(:,11); % distance in pixels between LEDs
datai = [rx,ry,gx,gy,spks,angle,distance]; % create array with all variables

% find red LED non-detects
dataFiltx = find(datai(:,1) > xmin & datai(:,1) < xmax & datai(:,2) < ymax & datai(:,1) > ymin);
rFILT = datai(dataFiltx,:);

% find green LED non-detects
dataFiltxy = find(rFILT(:,3) > xmin & rFILT(:,3) < xmax & rFILT(:,4) < ymax & rFILT(:,4) > ymin);
rgFILT = rFILT(dataFiltxy,:);

% find Min and Max LED distance
dataFiltxLED = find(rgFILT(:,7) > minLED & rgFILT(:,7) < maxLED);
rgmmFILT = rgFILT(dataFiltxLED,:);

% extract x y coords and smooth 
rxs = runline(rgmmFILT(:,1),5,1); % smooth with 5pt window and 1pt step (from Chronux toolbox)
rys = runline(rgmmFILT(:,2),5,1);

% velocity of rat from smoothed xy data
vel_x = diff(rxs); % vel units are pixels/frame
vel_y = diff(rys);
vel_abs = sqrt(vel_x.^2 + vel_y.^2); % scalar length of velocity vector = "scalar velocity" in pixels/frame
vel_cmPerSec = vel_abs * 2.1 * sampleRate; % ~2.1 cm/pixel for Taube lab according to Shawn W  

% find xy coords and spikes when rat is above velocity threshold
iVel = find(vel_cmPerSec >= 0); % 5cm/sec is used by Stensola et al 2012
spksf = rgmmFILT(:,5);
dataVEL = [rxs rys spksf];
dataVEL = dataVEL(iVel,:);

% extract spike locations from smoothed and velocity filtered data
spksi = find(dataVEL(:,3) == 1); 
spks_VEL = dataVEL(spksi,:);

% create array for binning data
nBins = 26; % adjusted so the bins are ~5x5cm as in Koenig et al 2011
MinY = min(dataVEL(:,2));
MaxY = max(dataVEL(:,2));
MinX = min(dataVEL(:,1));
MaxX = max(dataVEL(:,1));
edges{1} = linspace(MinY, MaxY, nBins+1);
edges{2} = linspace(MinX, MaxX, nBins+1);

% bin occupancy data
occMatrix = [dataVEL(:,2),dataVEL(:,1)];
Omatrix = hist3(occMatrix,'Edges',edges);
Omatrix(1,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/sampleRate;

% bin spike data
spikeMatrix = [spks_VEL(:,2), spks_VEL(:,1)];
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
[SmoothRateMap1] = smooth(FilledRateMatrix,1,5,1,5); % smooth with 5x5pt window and 1pt step

% end

data = importdata(ReadData{2});
% extract coords, spikes, angle, and direction from ReadData output
rx = data.data(:,2); % red x-coords
ry = data.data(:,3); % red y-coords
gx = data.data(:,4); % green x-coords
gy = data.data(:,5); % green y-coords
spks = data.data(:,6); % spikes
angle = data.data(:,10); % angle in radians
distance = data.data(:,11); % distance in pixels between LEDs
datai = [rx,ry,gx,gy,spks,angle,distance]; % create array with all variables

% find red LED non-detects
dataFiltx = find(datai(:,1) > xmin & datai(:,1) < xmax & datai(:,2) < ymax & datai(:,1) > ymin);
rFILT = datai(dataFiltx,:);

% find green LED non-detects
dataFiltxy = find(rFILT(:,3) > xmin & rFILT(:,3) < xmax & rFILT(:,4) < ymax & rFILT(:,4) > ymin);
rgFILT = rFILT(dataFiltxy,:);

% find Min and Max LED distance
dataFiltxLED = find(rgFILT(:,7) > minLED & rgFILT(:,7) < maxLED);
rgmmFILT = rgFILT(dataFiltxLED,:);

% extract x y coords and smooth 
rxs = runline(rgmmFILT(:,1),5,1); % smooth with 5pt window and 1pt step (from Chronux toolbox)
rys = runline(rgmmFILT(:,2),5,1);

% velocity of rat from smoothed xy data
vel_x = diff(rxs); % vel units are pixels/frame
vel_y = diff(rys);
vel_abs = sqrt(vel_x.^2 + vel_y.^2); % scalar length of velocity vector = "scalar velocity" in pixels/frame
vel_cmPerSec = vel_abs * 2.1 * sampleRate; % ~2.1 cm/pixel for Taube lab according to Shawn W  

% find xy coords and spikes when rat is above velocity threshold
iVel = find(vel_cmPerSec >= 0); % 5cm/sec is used by Stensola et al 2012
spksf = rgmmFILT(:,5);
dataVEL = [rxs rys spksf];
dataVEL = dataVEL(iVel,:);

% extract spike locations from smoothed and velocity filtered data
spksi = find(dataVEL(:,3) == 1); 
spks_VEL = dataVEL(spksi,:);

% create array for binning data
nBins = 26; % adjusted so the bins are ~5x5cm as in Koenig et al 2011
MinY = min(dataVEL(:,2));
MaxY = max(dataVEL(:,2));
MinX = min(dataVEL(:,1));
MaxX = max(dataVEL(:,1));
edges{1} = linspace(MinY, MaxY, nBins+1);
edges{2} = linspace(MinX, MaxX, nBins+1);

% bin occupancy data
occMatrix = [dataVEL(:,2),dataVEL(:,1)];
Omatrix = hist3(occMatrix,'Edges',edges);
Omatrix(1,:) = [];
Omatrix(:,end) = [];
occ = Omatrix/sampleRate;

% bin spike data
spikeMatrix = [spks_VEL(:,2), spks_VEL(:,1)];
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
[SmoothRateMap2] = smooth(FilledRateMatrix,1,5,1,5); % smooth with 5x5pt window and 1pt step

SpatialCrossCorr = corr2(SmoothRateMap1,SmoothRateMap2);

% plot filled and smoothed rate map 
subplot(2,1,1), h = pcolor(SmoothRateMap1);
axis square tight
hold on
colorbar;
box off
set(h, 'EdgeColor', 'none');
hold on
subplot(2,1,2), h = pcolor(SmoothRateMap2);
axis square tight
hold on
colorbar;
box off
set(h, 'EdgeColor', 'none');

keep SpatialCrossCorr
