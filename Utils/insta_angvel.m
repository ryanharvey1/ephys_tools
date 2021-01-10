function anglevel = insta_angvel(theta,samplerate,smooth_factor)
% insta_angvel: calc instant angular velocity
%
% Input:
%           theta: angles in radians
%           samplerate: samples per second
%           smooth_factor: integer value for which to divide sample rate
%           (i.e. samplerate/smooth_factor). 
%
% Output:
%           anglevel: angular velocity in deg/sec
%
% Ryan H, Laura B; Revised by LB 01/21
%

angdist = circ_dist(theta(1:end-1),theta(2:end)); % -pi:pi
anglevel = rad2deg(angdist*samplerate);

% apply Savitzky-Golay filtering to remove artifacts due to NaNs or
% Tracking flips
frames_to_smooth = nearest_odd(samplerate/smooth_factor);
anglevel = sgolayfilt(anglevel,3,frames_to_smooth); 

end

% to set the number of frames to smooth - must be odd for Savitzky-Golay
function y = nearest_odd(x)
y = 2*floor(x/2)+1;
end


