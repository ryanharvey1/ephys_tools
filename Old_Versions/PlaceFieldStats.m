% creates summary of data for output from Ryan's data
% 
% Input: 
%       - path to raw data output from Ryan Yoder's programs
%
% Output:
%       - rate map
%       - rate map stats (sparsity, information content)
% 
% Created by Ben C June 2015; Updated August 2016

addfolders = 'F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\Cell analysis';
addpath(genpath(addfolders)); % added for FindFiles

% load rate map array
% path = '/Users/benclark/Dropbox/Place Cell Analysis/Raw Data/RH026 10-28-2013/';
path = 'D:\Place_Cell_Data\Place Cells - Tilted Mice\2013-08-01_08-56-12 - SR012 - cells RY analyzed\Ascii New';
% SpikeFile = FindFiles('TT2- RH_SS_09.txtrmap1.txt', 'StartingDirectory', path);
SpikeFile = FindFiles('TT1-RY 2_SS_01.txtrmap1.txt', 'StartingDirectory', path);
data = importdata(SpikeFile{1});

% extract coords, spikes, angle, and direction from ReadData output
RateMap = data.data; % rate map
x_up = (numel(RateMap(1,:)));
x_low = x_up-27;
y_up = (numel(RateMap(:,1)));
y_low = y_up-27;
RateMap(1:y_low, :) = [];
RateMap(:,1:x_low) = [];

% calculate rate map statistics
xbins = numel(RateMap(1,:));
ybins = numel(RateMap(:,1));
r = reshape(RateMap,xbins*ybins,1); % reshape data into column
r(isnan(r))=[];
Peak_Rate = max(r);
r_25perc = Peak_Rate*0.25;
r_25percAll = find(r > r_25perc);
NumbActiveBins = numel(r_25percAll);
AllBins = numel(r);
PercActiveBins = (NumbActiveBins/AllBins)*100;

[y_max, x_max] = find(RateMap == Peak_Rate);
x_top = 27 - x_max;
x_bottom = x_max - 1;
y_top = 27 - y_max;
y_bottom = y_max - 1;
distanceAll = [x_top x_bottom y_top y_bottom];
distanceMin = min(distanceAll);
distance_to_wall = distanceMin*2.25;

% plot raw rate map 
figure(1), h = pcolor(RateMap);
colormap(jet);
axis square tight
hold on
set(h, 'EdgeColor', 'none');
colorbar;
box off
set(gca,'YDir','reverse');

keep RateMap PercActiveBins distance_to_wall Peak_Rate