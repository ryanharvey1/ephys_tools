function [ E ] = Eccentricity( X,Y )
%Finds the eccentricity (of ellipse) in the Morris Water Maze
%   X = all x-coordinates of the trajectory
%   Y = all y-coordinates of the trajectory
[a,b,~]=min_encl_ellipsoid(X,Y);
E=sqrt(1-(b^2/a^2));
end

