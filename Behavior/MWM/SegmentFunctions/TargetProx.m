function [ TP ] = TargetProx( X,Y,Px,Py,r )
%Function finds the percentage of a trajectory that is within a radius of
%the platform on the Morris water maze
%   - r is the radius of the platform
%   - X is the x coordinates of every point along the trajectory
%   - Y is the y coordinates of every point along the trajectory
%   - Px is the center of the platform's x coordinate
%   - Py is the center of the platform's y coordinate
TP=(numel(find((sqrt((Px - X).^2+(Py - Y).^2))<(6*r))))/numel(X);

end

