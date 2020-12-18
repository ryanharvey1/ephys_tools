function anglevel = insta_angvel(theta,samplerate)
% insta_angvel: calc instant angular velocity
%
% Input:
%           theta: angles in degrees
%           samplerate: sample rate
%
% Output:
%           anglevel: angular velocity in angles/sec
%
% Ryan H, Laura B; updated 10/26/2020
%

angdist = circ_dist(theta(1:end-1),theta(2:end));
anglevel = angdist*samplerate;

end




