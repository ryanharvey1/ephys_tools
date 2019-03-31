% TuningCurveCrossCorr
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
rxf = rgmmFILT(:,1);
ryf = rgmmFILT(:,2);
rxs = runline(rxf,5,1); % smooth with 10pt window and 1pt step (based on Valente et al 2007 PloS One)
rys = runline(ryf,5,1);

% extract spike locations
spksf = rgmmFILT(:,5);
spksi = find(spksf == 1); 
spks_xs = rxs(spksi);
spks_ys = rys(spksi);
spks_xf = rxf(spksi);
spks_yf = ryf(spksi);

% extract angle and convert to degrees
anglef = rgmmFILT(:,6);
degHeading = circ_rad2ang(anglef);
dSpks = degHeading(spksi);
aSpks = anglef(spksi);
timestamps = 1:length(degHeading);
timestamps = timestamps';
tSpks = timestamps(spksi);

% number of directional bins
dBins = 60;

% Binning of directional data can also be done using histc function
for i = 1:dBins;
    % number of instances where the head of the rat was found in each of the 60 possible orientations
    ListOrientation = find(anglef >= ((((2*pi)/dBins)/2)+(i-1)*2*pi/dBins) & anglef < ((((2*pi/dBins)/2)+(i)*2*pi/dBins)));
        if length(ListOrientation) < 1; % if the number is less than 1 (i.e., 0), it is hard set to 1 to avoid division by 0
           nOrientation(i) = 1; 
           nSpikesOrientation(i) = 0; % 0 is assined to the number of spikes for this orientation
        else
           nSpikesOrientation(i) = sum(spksf(ListOrientation));
           nOrientation(i) = length(ListOrientation);
    end
end

% transformed values from 1/60th of a sec to seconds
nOrientation2 = nOrientation./sampleRate; % 60ms for the time in s

% calculates the spikes/sec for each 6° bin
BinsNbSpikes = nSpikesOrientation./nOrientation2;

% trace the graph of the discharge rate / direction to the cells 1
BinsAngle = (0+((2*pi/(dBins)/2)):2*pi/(dBins):(2*pi)-((2*pi/(dBins)/2)));
BinsAngle3 = BinsAngle*dBins;


% 
% 
% SpatialCrossCorr = corr2(SmoothRateMap1,SmoothRateMap2);
% 
% % plot filled and smoothed rate map 
% subplot(2,1,1), h = pcolor(SmoothRateMap1);
% axis square tight
% hold on
% colorbar;
% box off
% set(h, 'EdgeColor', 'none');
% hold on
% subplot(2,1,2), h = pcolor(SmoothRateMap2);
% axis square tight
% hold on
% colorbar;
% box off
% set(h, 'EdgeColor', 'none');
% 
% keep SpatialCrossCorr
