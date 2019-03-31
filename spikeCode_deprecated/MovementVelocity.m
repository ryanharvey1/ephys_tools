% MovementVelocity creates summary of movement charactersitics
%
% Input: 
%       - raw data output from Labview
%
% Output:
%       - 
% 
% Created by Ben C June 2014

% identify path to data files
path = '/Users/bjclark/Desktop/Dropbox/EC_analysis/__RawData_&_Summary/__control_results_files_MEC';
ReadData = FindFiles('bc217*.txt', 'StartingDirectory', path);

for i = 1:length(ReadData);

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

% load xy data
data = importdata(ReadData{i});

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

meanLinearVel = mean(vel_cmPerSec);

% data and images to retain
keep('path', 'ReadData', 'sampleRate', 'xmin', 'xmax', 'ymin', 'ymax', 'minLED', 'maxLED', 'dBins', 'i', 'meanLinearVel');
[filepath, filename] = fileparts(ReadData{i});
save([filepath filesep filename '_Velocity.mat']);

end