function [xc,yc,r] = fitcircle(x,y)

% [xc,yc,r] = fitcircle(x,y)
% 
% INPUTS:
%     x, y -- x,y coordinates of a rat running on a circular track
%
% OUTPUTS:
%     xc, yc -- x and y coordinate of the center of the track
%     r -- radius of the track
%
% AT 2011

xx = x.*x;
yy = y.*y;
xy = x.*y;

A = [sum(x)  sum(y)  length(x);
     sum(xy) sum(yy) sum(y);
     sum(xx) sum(xy) sum(x)];
B = [-sum(xx+yy);
     -sum(xx.*y+yy.*y);
     -sum(xx.*x+xy.*y)];

a = A \ B;
xc = a(1) / -2;
yc = a(2) / -2;
r  = sqrt((a(1)^2 + a(2)^2) / 4 - a(3));