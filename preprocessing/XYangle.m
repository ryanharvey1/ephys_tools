function [ angle ] = XYangle(X,Y)
%XYangle Obtains angle of heading from XY path*
%  Input
%           X: x dim of data in a column vector
%           Y: y dim of data in a column vector
%  Output
%       angle: single column vector of angle of heading (will be -1 the length of input vector)
% 
% *Main calculation: angle= atan2d(y2-y1,x2-x1) + 360*((y2-y1)<0);
% 
% Ryan Harvey 6/2/17
% 
for i=1:length(X)
    if i+1>length(X);break;end
    angle(i,1)= atan2d(Y(i+1,1)-Y(i,1),X(i+1,1)-X(i,1)) + 360*((Y(i+1,1)-Y(i,1))<0);
end


