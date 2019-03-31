function [ vel_cmPerSec,vel_abs,pixelDist ] = InstaVel(Position,linear_track,track_length,sampleRate)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% velocity of rat from smoothed xy data
% scalar length of velocity vector = "scalar velocity" in pixels/frame
if isequal(linear_track,'yes')
    vel_abs = abs(diff(Position(:,1)));
else
    vel_abs = sqrt(diff(Position(:,1)).^2 + diff(Position(:,2)).^2);
end
pixelDist = track_length/(range(Position(:,1)));
vel_cmPerSec = vel_abs * pixelDist * sampleRate;

end

