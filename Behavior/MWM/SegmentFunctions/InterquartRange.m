function [ IQR ] = InterquartRange( X,Y,Mx,My,R )
%Function finds the interquartile range of the distance to the center of a
%trajectory in the morris water maze
%   - X is the x coordinates of each point on the trajectory
%   - Y is the y coordinates of each point on the trajectory
%   - Mx is the x coordinate of the center of the MWM
%   - My is the y coordinate of the center of the MWM
%   - R is the radius of the MWM
IQR=(quantile(sqrt((Mx-X).^2+(My-Y).^2),.75)-quantile(sqrt((Mx-X).^2+(My-Y).^2),.25))/R;

end

