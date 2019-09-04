function [ TL ] = TrajLeng( X,Y )
%Finds the length of the trajectory of a rat in the Morris Water Maze
%  Inputs
%   X=the x coordinates of the trajectory
%   Y=the y coordinates of the trajectory
 
TL=sum(sqrt(diff(X).^2+diff(Y).^2));
end

