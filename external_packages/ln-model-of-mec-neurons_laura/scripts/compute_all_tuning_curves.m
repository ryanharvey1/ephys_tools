%% Description
% This will compute the firing rate tuning curves for position, head
% direction, running speed, and theta phase.
%
% Updated to include egocentric bearing and distance from maze boundary.
% Laura Berkowitz 02/2020

% take out times when the animal ran >= 50 cm/s
posx_c(too_fast) = []; posy_c(too_fast) = []; 
direction(too_fast) = [];
% speed(too_fast) = [];
phase(too_fast) = [];
ego_bearing(too_fast) = [];
% distances(too_fast) = [];

% compute tuning curves for position, head direction, speed, and theta phase
[pos_curve] = compute_2d_tuning_curve(posx_c,posy_c,smooth_fr,n_pos_bins,0,boxSize);
[hd_curve] = compute_1d_tuning_curve(direction,smooth_fr,n_dir_bins,0,2*pi);

% [speed_curve] = compute_1d_tuning_curve(speed,smooth_fr,n_speed_bins,0,50);
[theta_curve] = compute_1d_tuning_curve(phase,fr,n_theta_bins,0,2*pi);

[ego_curve] = compute_1d_tuning_curve(ego_bearing,smooth_fr,n_ego_bins,0,360);
% [dist_curve] = compute_1d_tuning_curve(distances,smooth_fr,n_ego_bins,0,dist_range(1,2)*2);