
%%%%%%%% DATA EXTRACTION & PREPROCESSING from Neuralynx Files %%%%%%%%%%%%%

% Ryan H and Ben C May 2016

% load path / video file / spike files
path ='D:\Place_Cell_Data\PAE_Rat\RH16\2016-07-02_18-46-42'; % For .nvt file
% add needed functions to path 
addfolders = 'F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\spikeCode';
addpath(genpath(addfolders));  
addfolders = 'F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\Cell analysis';
addpath(genpath(addfolders));  
addfolders = 'F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\chronux_2_11';
addpath(genpath(addfolders));  

VTfile = FindFiles('*.nvt', 'StartingDirectory', path); 
try
    SpikeFile = FindFiles('*SS*.txt', 'StartingDirectory', path); % For Spikesort 3D
catch
end
if exist('SpikeFile')==2; % For Spikesort 3D
    SS3D=2;
end
if  exist('SpikeFile')==1
    NSMA_Toolbox = ('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox'); % add NSMA to path
    addpath(genpath(NSMA_Toolbox));
    mclustpath=cd([path,'p','\TT']); % sets path to p file & looks within TT folder
    ls '*.t' % disp .t files
    tfile=FindFiles('*.t');

for i=1:length(tfile); % Counts # of .t files (# of clusters) and sets iteration length
    SS3D=1; % For MClust
    NSMA_Toolbox = ('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\NSMA_Toolbox'); % add NSMA to path
    addpath(genpath(NSMA_Toolbox));
    S =LoadSpikes(tfile); % loads .t into ts cells
    t = Data(S{i}); % converts ts cells to numerals 
    SpikeFile=100*t; % conversion from MClust output to our current code (subject to change)
    %SpikeFile = num2cell(SpikeFile); % converts non-cell array object to cell array
    rmpath(genpath(NSMA_Toolbox)); % Remove NSMA from path (different Nlynx code)

% for i = 1:length(SpikeFile);  % For Spikesort 3D (commented out, still working on it.)

% video frame rate in Hz
sampleRate = 30;

% extract video data with mex
% Nlx2MatVT_v3 for mac compiler *** Nlx2MatVT for PC compiler 
[VTTimeStamps, ExtractedX, ExtractedY, ExtractedAngle] = Nlx2MatVT(VTfile{1}, [1 1 1 1 0 0], 0, 1);                                                      
ts_video = VTTimeStamps';
angle_video = ExtractedAngle';
data_video = [ts_video, ExtractedX', ExtractedY', angle_video];

%Filter out tracking issues
data_video1 = (data_video(:,2)>10); % can be changed based on maze config
data_video = data_video(data_video1,:);
data_video1 = (data_video(:,3)>10); % can be changed based on maze config
data_video = data_video(data_video1,:);

% Smoothing data
x_smooth = runline(data_video(:,2),10,1); 
y_smooth = runline(data_video(:,3),10,1);
data_video_smoothfilt = [data_video(:,1) x_smooth y_smooth data_video(:,4)];

%Interpolate 
if SS3D==2
    s1 = importdata(SpikeFile{i});       % for spikesort 3D text files
end
if SS3D==1
    s1 = SpikeFile;                    % for MClust .t files
end
ts_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,1), s1, 'nearest');
% x_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,2), s1, 'nearest');
% y_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,3), s1, 'nearest');
% angle_spike = interp1(data_video_smoothfilt(:,1), data_video_smoothfilt(:,4), s1, 'nearest');

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
vel_cmPerSec = vel_abs * 2.1 * sampleRate; % ~2.1 cm/pixel for Taube lab according to Shawn W  

% find xy coords and spikes when rat is above velocity threshold
iVel = find(vel_cmPerSec >= 0); % 5cm/sec is used by Stensola et al 2012
data_video_smoothvelfilt = data_video_smoothfilt(iVel,:);

% extract spike locations from smoothed and velocity filtered data
spksi = find(data_video_smoothvelfilt(:,5) == 1); 
spks_VEL = data_video_smoothvelfilt(spksi,:);

% create array for binning data
nBinsx = 20; % adjusted so the bins are ~5x5cm as in Koenig et al 2011
nBinsy = 8;
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

% NEW find pixels with zero occupancy- these will get blanked out out in the
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

% NEW
SmoothRateMap(Occ1ZeroLogical) = NaN;
FilledRateMatrix(Occ1ZeroLogical) = NaN;

% plot filtered and smoothed rate map 
figure (i), subplot(2,1,1), plot(data_video_smoothvelfilt(:,2), data_video_smoothvelfilt(:,3), 'LineWidth', 2, 'color', [0,0,0]+0.8);
hold on;
scatter(spks_VEL(:,2), spks_VEL(:,3), 45, 'filled', 'k');
box off
axis off
title('Spike (black dots) on Path (gray)');

% plot filled and smoothed rate map 
% NaNsForSmoothedR = NaN(1,nBins+1);
% NaNsForSmoothedC = NaN(nBins,1);
% SmoothRateMapGraph = [SmoothRateMap,NaNsForSmoothedC;NaNsForSmoothedR];
figure (i);, subplot(2,1,2), h = pcolor(SmoothRateMap);
colormap jet
axis off
hold on
% colorbar
box off
set(h, 'EdgeColor', 'none');
% set(gca,'YDir','reverse');
title('Smoothed Rate Map');

% images / data to retain
if SS3D==2
%     cd(path);
    [filepath, filename] = fileparts(SpikeFile{i});
    save([filepath filesep filename '_pathProperties.mat']);
    print(figure (i), '-djpeg', filename);
    close(figure (i)); 
end
if SS3D==1
    try % temporary while testing
    cd(mclustpath);
    [filepath, filename] = fileparts(tfile{i});
    save([filepath filesep filename '_pathProperties.mat']);
    cd(mclustpath);
    print(figure (i), '-djpeg', filename);
    close(figure (i)); 
    catch
    end
    keep('path','VTfile','SpikeFile','mclustpath','tfile')
end
end
end
% end
% end




