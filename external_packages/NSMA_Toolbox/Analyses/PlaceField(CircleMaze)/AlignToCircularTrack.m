function [V,W] = AlignToCircularTrack(X,Y,xc,yc)

% AlignToCircularTrack  Reparameterizes x,y position data to degrees on circle
%
% [V,W] = AlignToCircularTrack(X,Y,xc,yc)
%
% INPUTS: 
%   X,Y - (c)tsd or straight arrays of x,y coordinates
%   xc,yc - x,y coordinates of center of circular track
% OUTPUTS:
%   V,W - (c)tsd or straight arrays of data, depending on form of inputs:
%       V = position along circle (degrees), W = radial distance from center 
%
% Called within ParameterizeToCircularTrack function
%
% ADR 1998, last modified '04 by MN


x = Data(X);
y = Data(Y);

%--------------------
% Project
%--------------------
W0 = sqrt((x - xc).* (x - xc) + (y - yc) .* (y - yc));
V0 = 180/pi * (atan2(y - yc, x - xc));

%--------------------
% restrict 0 <= V0 < 360
%--------------------
f = find(V0 < 0);
V0(f) = V0(f) + 360;
f = find(V0 >= 360);
V0(f) = V0(f) - 360;

%--------------------
% Build standard data structure
%--------------------
V = tsd(Range(X,'ts'), V0);
W = tsd(Range(X,'ts'), W0);

