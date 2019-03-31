function [V,W] = ParamToHorseshoeTrack(X, Y)

% [V, W] = ParamToHorseshoeTrack(X,Y)
% Parameterize X,Y position data of a rat running on a circular track in a
% horseshoe pattern to V (angle around circle) and W (radial distance from center) coordinates
%
% INPUTS:
%     X, Y -- ctsd or tsd or straight arrays of data
%
% OUTPUTS:
%     V, W -- ctsd or tsd or straight arrays of data
%     V = position along circle, W = radial distance from center
%     output is of same type as input
%
% AT 2011: fits a circle to the track to automatically find the center of
%  the track coordinates, and then finds the turn around point and shifts
%  the V-coordinate to it.
% calls AlignToCircularTrack and fitcircle

%------------------------
x = Data(X);
y = Data(Y);

% fit a circle to the data to find the centre of the track
[xc,yc,r] = fitcircle(x, y);

% scale the first data point to the track and use it as a starting point
p1 = atan2(y(1) - yc, x(1) - xc);
x0 = xc + r*cos(p1);
y0 = yc + r*sin(p1);

% test code to display output
t = linspace(0,2*pi,360);
x1 = xc + r*cos(t);
y1 = yc + r*sin(t);
figure; hold; xlabel('X'); ylabel('Y'); axis equal;
plot(x, y, 'k-');              % plot the data points
plot(x1,y1,'r','LineWidth',2); % plot the fitted circle
plot(x0,y0,'b*');              % plot the starting point

[V,W] = AlignToCircularTrack(X,Y, xc, yc, x0, y0);

% align the parameterized data to the start of the track and centre it
V_aligned = unwrap(Data(V));
V_floor = (max(V_aligned) + min(V_aligned) - 360)/2;
V_aligned = V_aligned - V_floor;

% rebuild V with the aligned data
switch class(X)
   case 'tsd'
      V = tsd(Range(V,'ts'), V_aligned);
   case 'ctsd'
      V = ctsd(StartTime(V), DT(V), V_aligned);
   otherwise
      V = V_aligned;
end
