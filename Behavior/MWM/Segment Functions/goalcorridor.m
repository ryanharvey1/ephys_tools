function [ GC ] = goalcorridor ( X,Y,Px,Py )
%Finds the percentage of points of a rat's trajectory within the goal corridor in the
%Morris Water Maze
%Goal corridor is between two lines, one 20 degrees above and the other below the straight
%line from the starting poing to the platform
%   - X is the x coordinates of every point along the trajectory
%   - Y is the y coordinates of every point along the trajectory
%   - Px is the center of the platform's x coordinate
%   - Py is the center of the platform's y coordinate
%finding the degree of the straight line from the starting point to the
%platform
Sx=X(1,1);
Sy=Y(1,1);
heading=XYangle([Sx; Px],[Sy; Py]);   
alt=sqrt((Sx-Px)^2+(Sy-Py)^2); 
lineLength = alt+10;
angle1 = heading+20;
angle2 = heading-20;
x(1) = Sx;
y(1) = Sy;
x(2) = x(1) + lineLength * cosd(angle1);
y(2) = y(1) + lineLength * sind(angle1);

x(3) = x(1) + lineLength * cosd(angle2);
y(3) = y(1) + lineLength * sind(angle2);
inGC=inpolygon(X,Y,x',y');
GC=(sum(inGC)/size(X,1));

end
