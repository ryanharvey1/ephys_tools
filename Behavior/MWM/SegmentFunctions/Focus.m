function [ F ] = Focus( X,Y )
%Finds the focus of a rat's trajectory of a Morris Water Maze
%   X = x-coordinates of the trajectory points
%   Y = y-coordinates of the trajectory points
[a,b,~]=min_encl_ellipsoid(X,Y);
F=1-4*a*b/(TrajLeng(X,Y)^2);
end

