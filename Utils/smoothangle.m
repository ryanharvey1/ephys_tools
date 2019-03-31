function [thetaout]=smoothangle(thetain,window)
% smoothangle
% Input:
%       thetain: angle in radians 
%       window: smoothing window
% Ouput:
%       thetaout: smoothed angles in radians (-pi to pi)
% 
% Ryan E Harvey 2018


% take difference in discontinuous angles in radius
delta_theta=diff(thetain); 

% derivative of circular data
delta_theta=atan2(sin(delta_theta), cos(delta_theta)); 

% continuized data
thetaout= [thetain(1);cumsum(delta_theta) + thetain(1)];

% smooth with gaussian at a moving window 
thetaout=smoothdata(thetaout,'gaussian',window); 

% bring back to radians bound by [-pi, pi)
thetaout=atan2(sin(thetaout),cos(thetaout)); 
end



