function [ MD2 ] = MedD2center( X,Y,Mx,My,r )
%Function finds the median distance to the center of a trajectory in the
%Moriss Water maze
%   - X is the x coordinates of each point in the trajectory
%   - Y is the y coordinates of each point in the trajectory
%   - Mx is the x coordinate of the center of the MWM
%   - My is the Y coordinate of the center of the MWM
%   - r is the radius of the MWM
MD2 = median(sqrt((Mx-X).^2+(My-Y).^2))/r;


end

