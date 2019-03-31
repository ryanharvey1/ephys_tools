function [ z ] = armtoradians(y)

% turns arms of the 8 arm radial arm maze into radians
% for i=size(y);
y(y==1) = 0;
y(y==2) = 45;
y(y==3) = 90;
y(y==4) = 135;
y(y==5) = 180;
y(y==6) = 225;
y(y==7) = 270;
y(y==8) = 315;
[z]=circ_ang2rad(y');
% end
end

