% oneCell_SpikesRateMap_Ryan creates summary of data for output from Ryan's
% data
% 
% Input: 
%       - raw data output from Ryan Yoder's programs
%
% Output:
%       - rate map
%       - rate map stats (sparsity, information content)
% 
% Created by Ben C Jan 2015

% identify path to data files
path = '/Users/bjclark/Desktop/Dropbox/place cell analysis/Raw Data/RH048/9-11-2013';
ReadData = FindFiles('TT1 - SK_SS_02.txtrmap1.txt', 'StartingDirectory', path);

% load rate map array
data = importdata(ReadData{1});

% extract coords, spikes, angle, and direction from ReadData output
RateMap = data.data; % rate map

% calculate rate map statistics
xbins = numel(RateMap(1,:));
ybins = numel(RateMap(:,1));
r = reshape(RateMap,xbins*ybins,1); % reshape data into column
r(isnan(r))=[];
r_max = max(r);
r_20perc = r_max*0.20;
r_20percAll = find(r > r_20perc);
NumbActiveBins = numel(r_20percAll);
AllBins = numel(r);
PercActiveBins = (NumbActiveBins/AllBins)*100;
% [Sparsity] = Sparsity(r'); % sparsity from NSMA toolbox
[x_max, y_max] = find(RateMap == r_max);
row_max = RateMap(x_max,:);
col_max = RateMap(:,y_max);
% dirBins = 0:6:180;
% ai = imrotate(RateMap,6,'nearest','crop');

% plot raw rate map 
figure(1), h = pcolor(RateMap);
axis square tight
hold on
colorbar;
box off
set(h, 'EdgeColor', 'none');
a = r > 0;
a = min(a);
[C, ~] = contour(RateMap,[a r_max],'r','LineWidth',2);
% set(gca, 'xdir','reverse');
% set(gca, 'ydir','reverse');

keep PercActiveBins C

% % data and images to retain
% [filepath, filename] = fileparts(ReadData{1});
% print(figure (1), '-djpeg', filename);