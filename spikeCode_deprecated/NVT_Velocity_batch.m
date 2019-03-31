%% load video file
path = '/Users/bjclark/Desktop/Dropbox/ExampleNlxRecording/2010-08-28_14-39-57';
VTfile = FindFiles('*.nvt', 'StartingDirectory', path);

for i = 1:length(VTfile);

% extract video data with mex
[VTTimeStamps, ExtractedX, ExtractedY] = Nlx2MatVT(VTfile{i}, [1 1 1 0 0 0], 0, 1);
data=[ExtractedX',ExtractedY'];

%% extract x y coords and smooth 
xs = runline(data(:,1),5,1); % smooth with 5pt window and 1pt step (from Chronux toolbox)
ys = runline(data(:,2),5,1);

%% velocity of rat from smoothed xy data
vel_x = diff(xs); % vel units are pixels/frame
vel_y = diff(ys);
xCM = 120/minus(max(data(:,1)),min(data(:,1))); % Calculate number of pixels per cm
vel_xCM = vel_x.*xCM; % convert from pixels to cm between coordinates
yCM = 120/minus(max(data(:,2)),min(data(:,2))); % Calculate number of pixels per cm
vel_yCM = vel_y.*yCM; % convert from pixels to cm between coordinates
vel_abs = sqrt(vel_xCM.^2 + vel_yCM.^2); % calculate distance between two successive points
vel_cmPerSec = vel_abs/0.0166; % calculate speed in cm/sec (16.6 msec per bin)

%% calculate mean velocity per sliding window
Velocity=[];
win_width = 150;  %Sliding window width (sampling rate X sec)
slide_incr = 30;  %Slide for each iteration (sampling rate)
for j = 1:slide_incr:36000
    iVel = mean(vel_cmPerSec(j:j+win_width));  %Calculation for each window
    Velocity=[Velocity;iVel];
end

%% save data
[filepath, filename] = fileparts(VTfile{i});
save([filepath filesep filename '_velocity.mat']);
end
clear all