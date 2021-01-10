function [speed_grid,speedVec,velx] = ang_speed_map(angles,fs,nbins)
% Input: 
%   - angles: vector of head angles in radians
%   - fs: video sampling rate
%   - nbins: number of bins to bin over; 
% Output: 
%   - speed_grid: binary speed matrix (length of angles by number of speed bins).
%   - speedVec: edges for speed indicies              
%   - speed: angular velocity vector 
%
% Dependencies: 
%   - insta_angvel.m (ephys_tools>utils)

%compute velocity
velx = insta_angvel([angles(1);angles],fs,3); % smooth over 1/3 second
maxSpeed = 1000; velx(velx>maxSpeed) = maxSpeed; %send everything over 1000 deg/s to 1000 deg/s

speedVec = maxSpeed/nbins/2:maxSpeed/nbins:maxSpeed-maxSpeed/nbins/2;
speed_grid = zeros(numel(angles),numel(speedVec));

for i = 1:numel(angles)

    % figure out the speed index
    [~, idx] = min(abs(velx(i)-speedVec));
    speed_grid(i,idx) = 1;
    
end

return