function [x,y,in]=restrictMovement(x,y,direction)
% restrictMovement allows you to draw a line around xy coordinates in order
% to eliminate certain points you don't want...ie when the rat jumps out of
% maze or tracker errors. 
%
% Also, I have included a direction input argument so you have either
% restrict outside or inside points. This is valuable if you have a maze
% like a circular track where the rat could jump out away or towards the
% center of the maze. 
%
% Input         x,y: coordinates 
%         direction: 1 (default) to remove outside points; 0 to remove inside points
%         
%
% Output        x,y: retained coordinates
%                in: logical of which coordinates were retained (so you can index ts etc.)
%       
%
% Ryan Harvey

% check inputs
if nargin<3
    direction=1;
end

% set up figure
fig=figure;plot(x,y,'Color',[1,1,1,0.2]);hold on
title('Click around the points you want to keep')
xlabel('X')
ylabel('Y')
axis equal
darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
disp('PRESS "ENTER" TO EXIT')
i=1;

% let the user click around the coordinates
while true
    [X,Y]=ginput(1);
    if isempty(X)
        break
    end
    corners(i,:)=[X,Y];
    plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
    i=i+1;
end

% remove points outside or inside the shape
in=inpolygon(x,y,corners(:,1),corners(:,2));
if direction==1
    x=x(in);
    y=y(in);
else
    x=x(~in);
    y=y(~in);
    in=~in;
end
close
end