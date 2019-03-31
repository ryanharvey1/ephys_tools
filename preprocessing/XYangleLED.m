function [ angle ] = XYangleLED(X1,Y1,X2,Y2)
%XYangle Obtains angle of heading from XY path*
%  Input
%           X1: x coord of front LED (data in a column vector)
%           Y1: y coord of front LED (data in a column vector)
%           X2: x coord of back LED (data in a column vector)
%           Y2: y coord of back LED (data in a column vector)
%  Output
%       angle: single column vector of angle of heading
% 
% *Main calculation: angle= atan2d(y2-y1,x2-x1) + 360*((y2-y1)<0);
% 
% Ryan Harvey 6/2/17

    angle= deg2rad(atan2d(Y2-Y1,X2-X1) + 360*((Y2-Y1)<0));
end


