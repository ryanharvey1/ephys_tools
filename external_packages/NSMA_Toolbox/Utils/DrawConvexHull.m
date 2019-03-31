function [x,y] = DrawConvexHull

% DrawConvexHull  Allows user to draw convex hull on current axis; returns x,y points on hull
%
% [x,y] = DrawConvexHull
%
% INPUTS:
%       (none)
% OUTPUS:
%       x,y - lists of x and y coordinates of points on drawn hull
%
% ADR 1998, version V4.0, last modified '98 by ADR

% Status PROMOTED


[x,y] = ginput;
k = convhull(x,y);
x = x(k);
y = y(k);