% spikes_shuffle shifts the spike series relative to the time series
%
% Input: 
%       - raw data output from Labview
%
% Output:
%       - rSHUFF = mean vector length for shuffled spike series
%
% Created by Ben C 2014; Modified Ben C Feb 03 2017

% identify path to data files

clc, clear, close all

addpath(genpath('/Users/RyanHarvey/Dropbox/MATLAB/spikeCode'));
addpath(genpath('/Users/ryanharvey/Dropbox/MATLAB/CircStat2012a'));
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\spikeCode'));
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\CircStat2012a'));

path = 'F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\PAE_Project\ShuffleTestData';
ReadData = FindFiles('bc219*.txt', 'StartingDirectory', path);

% video frame rate in Hz
sampleRate = 60;

% maze parameters
xmin = 0;
xmax = 255;
ymin = 0;
ymax = 255;

% min and max LED distance parameters
minLED = 4;
maxLED = 39;

% number of directional bins
dBins = 60;
shuff_interv = dBins*20;

for i = 1:length(ReadData);
    
    data = importdata(ReadData{i}); % load xy data

    % extract coords, spikes, angle, and direction from ReadData output
    rx = data.data(:,2); % red x-coords
    ry = data.data(:,3); % red y-coords
    gx = data.data(:,4); % green x-coords
    gy = data.data(:,5); % green y-coords
    spks = data.data(:,6); % spikes
    angle = data.data(:,10); % angle
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

    % extract spike series and angle
    spksf = rgmmFILT(:,5);
    anglef = rgmmFILT(:,6);

    rSHUFF = [];
        for x = 1:400;
            spkSHUFFi = circshift(spksf,randi([-shuff_interv shuff_interv],1));
            for j = 1:dBins;
                ListOrientation = find(anglef >= ((((2*pi)/dBins)/2)+(j-1)*2*pi/dBins) & anglef < ((((2*pi/dBins)/2)+(j)*2*pi/dBins)));
                    if length(ListOrientation) < 1; 
                        nOrientation(j) = 1; 
                        nSpikesOrientation(j) = 0; % 0 is assined to the number of spikes for this orientation
                    else
                        nSpikesOrientation(j) = sum(spkSHUFFi(ListOrientation));
                        nOrientation(j) = length(ListOrientation); 
                    end
            end
            nOrientation2 = nOrientation./sampleRate;
            BinsNbSpikes = nSpikesOrientation./nOrientation2;
            BinsAngle = (0+((2*pi/(dBins)/2)):2*pi/(dBins):(2*pi)-((2*pi/(dBins)/2)));
            spksi = find(spkSHUFFi == 1);
            aSpks = anglef(spksi);
            r = circ_r(BinsAngle',BinsNbSpikes',circ_ang2rad(6));
            rSHUFF = [rSHUFF; r];            
        end
        
%         clear all non essential variables
%         keep('rSHUFF', 'corrSHUFF', 'path', 'ReadData', 'sampleRate', 'xmin', 'xmax', 'ymin', 'ymax', 'minLED', 'maxLED', 'dBins', 'shuff_interv', 'i');
        [filepath, filename] = fileparts(ReadData{i});
        save([filepath filesep filename '_SHUFF.mat']);  
end






