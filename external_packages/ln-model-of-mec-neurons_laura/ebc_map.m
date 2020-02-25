function [ego_grid,dist_grid, ego_bearing, distances, ego_bins,dist_bins, dist_range] = ebc_map(pos,hd, nbins)
% By Laura Berkowitz and Ryan Harvey February 2020
% input: 
%       pos: [x coords, y coords] in cm
%       hd: Head direction vector in degrees 
%       nbins: binsize 
%       boxSize: length of maze (assumes symmetrical arena)
%       mazeType: cylinder or box 
% output: 
%       ego_grid: length of pos by nBins. Identity matrix for egocentric
%       bearing. 
%       dist_grid: length of pos by nBins. Identify matrix for distance from boundary. 
%       ego_bearing: egocentric bearing across time. 
%       distances: distance from boundary across time. 
%       ego_bins: bin edges for egocentric bearing. 
%       dist_bins: bin edges for distance. 

% Dependencies: 
%   min_enclosed_elipse (path: ephys_tools\Utils)
%   fixNLXangle (path: ephys_tools\Utils)

% TO-DO (L.Berkowitz): 
%   add automated way to distinguish mazetype (e.g. cylinder or box) for boundary
%   calculations. 

hd = rad2deg(hd);
%find edges of environment
temp_pos = pos;
temp_pos(isnan(pos(:,1)),:) = [];
k = convhull(temp_pos);
xbound=temp_pos(k,1);
ybound=temp_pos(k,2);

clear temp_pos

% Obtain maximum radius
[ max_radius,~,center] = min_encl_ellipsoid(xbound,ybound);

dist_range = [0 max_radius/2];
% Create edges given maximum radius to create evenly spaced reference points
xbound=((cos(linspace(-pi,pi,360))*(max_radius+1)) + center(1,1))'; %add 1cm to max_radius to account for rearing
ybound=((sin(linspace(-pi,pi,360))*(max_radius+1))+ center(2,1))';

% calculate closest distance to boundary for all points
distances = zeros(length(pos),1); 
ego_ref = zeros(length(pos),1);

for i=1:length(pos)
    [ distances(i,1), ego_ref(i,1)]=min(sqrt(sum(bsxfun(@minus, [xbound ybound], [pos(i,1),pos(i,2)]).^2,2)));
end

% Calculate Egocentric Bearing
ego_bearing = wrapTo360(fixNLXangle(rad2deg(atan2(ybound(ego_ref,1) - pos(:,2),...
    xbound(ego_ref,1) - pos(:,1))),round(0.1667*30)));

ego_bearing=wrapTo360(hd - ego_bearing');

% define bin sizes for distance and egocentric heading
ego_bins = rad2deg(2*pi/nbins/2:2*pi/nbins:2*pi-2*pi/nbins/2);
dist_bins = (max_radius+5)/nbins/2:(max_radius+5)/nbins:(max_radius+5)-(max_radius+5)/nbins/2;

% initialize grid for distance and egocentric heading
ego_grid = zeros(length(pos),nbins);
dist_grid = zeros(length(pos),nbins);

% loop over positions to create identity matrices
for i = 1:length(ego_bearing)
    [~, idx] = min(abs(ego_bearing(i,1) - ego_bins));
    ego_grid(i,idx) = 1;
end

for i = 1:length(distances)
    [~, idx] = min(abs(distances(i,1) - dist_bins));
    dist_grid(i,idx) = 1;
end

end
