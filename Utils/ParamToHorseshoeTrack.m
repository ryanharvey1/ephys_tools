function [V,W] = ParamToHorseshoeTrack(x, y)

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
%
% Adapted from NSMA toolbox
% Ryan Harvey

% for code testing set to 1
figures=0;

% fit a circle to the data to find the centre of the track
% [xc,yc,r] = fitcircle(x, y);
[xc,yc,r]=betterfitcircle(x,y,figures);

% TO COMPENSATE FOR TRACKER ERRORS AND WHEN THE RAT LOOKS OVER THE BARRIER
% BETWEEN THE HORSESHOE, WE CAN SET THE STARTING POINT BASED ON WHERE THE
% XY COORDINATES ARE THE THINEST.
% THIS ASSUMES THE BARRIER IS AT THE Y MAX SECTION OF THE PLANE
xxs=unique(round(sort(x(~isnan(x)))));
xxs=xxs(round(length(xxs)/3):end-round((length(xxs)/3)));
in=zeros(length(xxs),1);
for i=1:length(xxs)
    [~,in(i,1)]= size(InterX([[xc';yc'],[xxs(i)';max(y)']],[x';y']));
end
[m,I]=nanmin(in);
P=InterX([[xc';yc'],[xxs(I)';max(y)']],[x';y']);
if m==0
    t = linspace(0,2*pi,360);
    x1 = xc + r*cos(t);
    y1 = yc + r*sin(t);
    P=InterX([[xc';yc'],[xxs(I)';max(y)*2']],[x1;y1]);
    y0=P(2,1);
else
    y0=median(P(2,:));
end
x0=xxs(I-1);


% test code to display output
if figures==1
    figure;
    t = linspace(0,2*pi,360);
    x1 = xc + r*cos(t);
    y1 = yc + r*sin(t);
    hold; xlabel('X'); ylabel('Y'); axis equal;
    plot(x, y, 'k-');              % plot the data points
    plot(x1,y1,'r','LineWidth',2); % plot the fitted circle
    plot(x0,y0,'b*');              % plot the starting point
    plot(xc,yc,'r*')               % plot center of circle
end

[V,W] = AlignToCircularTrack(x,y, xc, yc, x0, y0);

% align the parameterized data to the start of the track and centre it
V_aligned = unwrap(V);
V_floor = (max(V_aligned) + min(V_aligned) - 360)/2;
V_aligned = V_aligned - V_floor;

if max(V_aligned)<360
    V = V_aligned;
end

end

function [xCenter,yCenter,radius]=betterfitcircle(x,y,figures)
%  adapted from Image Analyst on matlab answers
% Image Analyst-
% "I solved it numerically by converting the data to a digital image,
% then using the Euclidean distance transform and finding the center of it.
% The max of the EDT is the center of the circle and its value there is
% the radius. If you need more than 3 digits of precision then you can 
% increase the size of the image from 1000 to 10000 or larger."

xtemp=x(~isnan(x));
ytemp=y(~isnan(y));

K = convhull(xtemp,ytemp);
xconvex=xtemp(K);
yconvex=ytemp(K);
if figures==1
    figure;
    fontSize = 12;
    subplot(2, 2, 1);
    plot(x,y,'.k')
    hold on
    plot(xconvex, yconvex, 'r', 'MarkerSize', 30);
    grid on;
    title('Original Points', 'FontSize', fontSize);
    axis image
end
% Enlarge figure to full screen.
% Make data into a 1000x1000 image.
xMin = min(xconvex);
xMax = max(xconvex);
yMin = min(yconvex);
yMax = max(yconvex);
scalingFactor = 1000 / min([xMax-xMin, yMax-yMin]);
x2 = (xconvex - xMin) * scalingFactor + 1;
y2 = (yconvex - yMin) * scalingFactor + 1;
mask = poly2mask(x2, y2, ceil(max(y2)), ceil(max(x2)));
% Display the image.
if figures==1
    p2 = subplot(2, 2, 2);
    imshow(mask);
    axis(p2, 'on', 'xy');
    title('Mask Image', 'FontSize', fontSize);
    axis image
end
% Compute the Euclidean Distance Transform
edtImage = bwdist(~mask);
% Display the image.
if figures==1
    p3 = subplot(2, 2, 3);
    imshow(edtImage, []);
    axis image
    axis(p3, 'on', 'xy');
end
% Find the max
radius = max(edtImage(:));
% Find the center
[yCenter, xCenter] = find(edtImage == radius);
% if multiple center points are located, find the mean between them. 
yCenter=mean(yCenter);
xCenter=mean(xCenter);

if figures==1
    % Display circles over edt image.
    viscircles(p3, [xCenter, yCenter], radius);
    % Display polygon over image also.
    hold on;
    plot(x2, y2, 'r', 'MarkerSize', 30, 'LineWidth', 2);
    title('Euclidean Distance Transform with Circle on it', 'FontSize', fontSize);
    % Display the plot again.
    subplot(2, 2, 4);
    plot(x,y,'.k')
    hold on
    plot(xconvex, yconvex, 'r', 'MarkerSize', 30);
    grid on;
    % Show the circle on it.
    hold on;
end
% Scale and shift the center back to the original coordinates.
xCenter = (xCenter - 1)/ scalingFactor + xMin;
yCenter = (yCenter - 1)/ scalingFactor + yMin;
radius = radius / scalingFactor;
if figures==1
    rectangle('Position',[xCenter-radius, yCenter-radius, 2*radius, 2*radius],'Curvature',[1,1],'EdgeColor','b');
    title('Original Points with Inscribed Circle', 'FontSize', fontSize);
    axis image
end

end
function [xc,yc,r] = fitcircle(x,y)

% [xc,yc,r] = fitcircle(x,y)
%
% INPUTS:
%     x, y -- x,y coordinates of a rat running on a circular track
%
% OUTPUTS:
%     xc, yc -- x and y coordinate of the center of the track
%     r -- radius of the track
%
% AT 2011

xx = x.*x;
yy = y.*y;
xy = x.*y;

A = [nansum(x)  nansum(y)  length(x);
    nansum(xy) nansum(yy) nansum(y);
    nansum(xx) nansum(xy) nansum(x)];
B = [-nansum(xx+yy);
    -nansum(xx.*y+yy.*y);
    -nansum(xx.*x+xy.*y)];

a = A \ B;
xc = a(1) / -2;
yc = a(2) / -2;
r  = sqrt((a(1)^2 + a(2)^2) / 4 - a(3));
end

function [V,W] = AlignToCircularTrack(X,Y,xc,yc,x0,y0, varargin)
%
% [V,W] = AlignToCircularTrack(X,Y,xc,yc,x0,y0, parameters)
%
% INPUTS:
%    X, Y -- ctsd or tsd or straight arrays of data
%    xc,yc -- center of track
%    x0,y0 -- point on track (used for orientation of coordinate system)
%
% OUTPUT:
%    V, W -- ctsd or tsd or straight arrays of data
%    V = position along circle, W = radial distance from center
%	output is of same type as input
%
% PARAMETERS
%    direction (default cw, also can be ccw)

% ADR 1998
%  version L4.1
%  status: PROMOTED

% v 4.1 added direction parameter

%--------------------
% Parms
direction = 'cw';
Extract_varargin;

%--------------------
% Check inputs
%--------------------

if (class(X) ~= class(Y))
    error('Both or none must be tsArrays.');
end
if (isa(X, 'tsd') | isa(X, 'ctsd'))
    CheckTS(X, Y);
    T0 = min(StartTime(X), StartTime(Y));
    T1 = max(EndTime(X), EndTime(Y));
    X = Restrict(X,T0,T1);
    x = Data(X);
    Y = Restrict(Y,T0,T1);
    y = Data(Y);
else
    x = X;
    y = Y;
end

%--------------------
% Get Orthonormal basis space
%--------------------

vx = atan2(y0 - yc, x0 - xc);

%--------------------
% Project
%--------------------
%calculates distance from each point to center of circle
W0 = sqrt((x - xc).* (x - xc) + (y - yc) .* (y - yc));
%calculates the angular difference between start and each point then converts to degrees
V0 = 180/pi * (atan2(y - yc, x - xc) - vx);

%--------------------
% restrict 0 < V0 < 360
%--------------------
f = find(V0 < 0);
V0(f) = V0(f) + 360;
f = find(V0 > 360);
V0(f) = V0(f) - 360;

if (strcmp(direction,'cw'))
    V0 = V0;
elseif (strcmp(direction, 'ccw'))
    V0 = 360 - V0;
else
    error('Unknown direction');
end

%--------------------
% Build standard data structure
%--------------------
switch class(X)                       % we know both X and Y are of the same class
    case 'tsd'
        V = tsd(Range(X,'ts'), V0);
        W = tsd(Range(X,'ts'), W0);
    case 'ctsd'
        V = ctsd(StartTime(X), DT(X), V0);
        W = ctsd(StartTime(X), DT(X), W0);
    otherwise
        V = V0;
        W = W0;
end
end
