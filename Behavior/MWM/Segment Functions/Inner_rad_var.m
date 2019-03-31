function [ IRS ] = Inner_rad_var( X,Y )
%Finds the Inner Radius Variation of a trajectory in the Morris Water Maze
%   X=all x-coordinates of the trajectory
%   Y=all y-coordinates of the trajectory
[~,~,c]=min_encl_ellipsoid(X,Y);
IRS=(quantile(sqrt((c(1)-X).^2+(c(2)-Y).^2),.75)-quantile(sqrt((c(1)-X).^2+(c(2)-Y).^2),.25))/mean(sqrt((c(1)-X).^2+(c(2)-Y).^2));

end

