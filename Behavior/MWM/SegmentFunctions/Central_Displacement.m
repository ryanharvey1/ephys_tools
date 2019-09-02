function [ CD ] = Central_Displacement( X,Y,Mx,My )
%Calculates the distance from the center of the minimum enclosing ellipse
%of a rat's trajectory in the Morris Water Maze to the center of the MWM.
%   X=the x coordinates of the trajectory
%   Y=the y coordinates of the trajectory
%   Mx=the MWM center x-coordinate
%   My=the MWM center y-coordinate
[~,~,c]=min_encl_ellipsoid(X,Y);
CD=sqrt((c(1)-Mx)^2+(c(2)-My)^2);
end

