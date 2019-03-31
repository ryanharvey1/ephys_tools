function [RateMap, PercActiveBins, distanceBins, Peak_Rate] = PlaceFieldStats(filename, x_low, y_low)
% creates summary of data for output from Ryan's data
% 
% Input: 
%       - path to raw data output from Ryan Yoder's programs
%
% Output:
%       - rate map
%       - rate map stats (sparsity, information content)
% 
% Created by Ben C June 2015

% load rate map array
data = importdata(filename);

% extract coords, spikes, angle, and direction from ReadData output
RateMap = data.data; % rate map
x_up = numel(RateMap(1,:));
y_up = numel(RateMap(:,1));

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
[y_max, x_max] = find(RateMap == r_max);
x_top = x_up - x_max;
x_bottom = x_max - x_low;
y_top = y_up - y_max;
y_bottom = y_max - y_low;

distanceAll = [x_top x_bottom y_top y_bottom];
distanceMin = min(distanceAll);
distanceBins = distanceMin;
Peak_Rate = r_max;

% plot raw rate map 
figure(1), h = pcolor(RateMap);
colormap(jet);
axis square tight
hold on
set(h, 'EdgeColor', 'none');
colorbar;
box off
scatter(x_max+0.5, y_max+0.5, 100, 'k', 's', 'LineWidth', 1.5);
set(gca,'YDir','reverse');
xlim([x_low x_up]);
ylim([y_low y_up]);
% set(gca,'XDir','reverse');

keep RateMap PercActiveBins distanceBins Peak_Rate