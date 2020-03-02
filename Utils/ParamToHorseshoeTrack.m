function [V,W] = ParamToHorseshoeTrack(x, y,varargin)

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

p = inputParser;
p.addParameter('figures',0); %for code testing
p.parse(varargin{:});
figures = p.Results.figures;


% fit a circle to the data to find the centre of the track
% [xc,yc,r] = fitcircle(x, y);
% [xc,yc,r]=betterfitcircle(x,y,figures);
% [xc,yc,r] = CircleFitByPratt([x,y]);
[xc,yc,r,~] = circfit(x,y);

% XY = [x,y];
% XY(any(isnan(XY),2),:) = []; % remove NaNs
% 
% [z, r, residual] = fitcircle(XY,'linear');
% xc = z(1);
% yc = z(2);

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
    subplot(4,2,[1,4])
    t = linspace(0,2*pi,360);
    x1 = xc + r*cos(t);
    y1 = yc + r*sin(t);
    hold; xlabel('X'); ylabel('Y'); axis equal;
    plot(x, y, 'k-');              % plot the data points
    plot(x1,y1,'r','LineWidth',2); % plot the fitted circle
    plot(x0,y0,'b*');              % plot the starting point
    plot(xc,yc,'r*')               % plot center of circle
    title('raw with circular fit')
end

[V,W] = AlignToCircularTrack(x,y, xc, yc, x0, y0);

% align the parameterized data to the start of the track and centre it
V_aligned = unwrap(V);
V_floor = (max(V_aligned) + min(V_aligned) - 360)/2;
V_aligned = V_aligned - V_floor;

if max(V_aligned)<360
    V = V_aligned;
end
if figures==1
    subplot(4,2,[5,6])
    plot(V,W)
    axis image
    title('unwrapped')
end

% detrend y axis
W = detrend_y(V,W);
if figures==1
    subplot(4,2,[7,8])
    plot(V,W)
    axis image
    title('detrend')
end
end

%% local functions below

function temp_w = detrend_y(V,W)
XY = [V, W];
XY(any(isnan(XY),2),:) = []; % remove NaNs

p = polyfit(XY(:,1),XY(:,2),8);

x1 = linspace(min(XY(:,1)),max(XY(:,1)),1000);
y1 = polyval(p,x1);

temp_w = W;
detrend_vec = (y1 - mean(y1));
for i = 1:1000-1
    temp_w(V > x1(i) & V < x1(i+1)) = temp_w(V > x1(i) & V < x1(i+1)) - detrend_vec(i);
end
end

function  [xc,yc,R,a] = circfit(x,y)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991, 

to_remove = isnan(x) | isnan(y);
x(to_remove) = [];
y(to_remove) = [];

x=x(:); y=y(:);
a=[x y ones(size(x))]\[-(x.^2+y.^2)];
xc = -.5*a(1);
yc = -.5*a(2);
R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end
function [xCenter,yCenter,radius] = CircleFitByPratt(XY)
%--------------------------------------------------------------------------
%  
%     Circle fit by Pratt
%      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
%      Computer Graphics, Vol. 21, pages 145-152 (1987)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this fit does not use built-in matrix functions (except "mean"),
%           so it can be easily programmed in any programming language
%
%--------------------------------------------------------------------------
XY(any(isnan(XY),2),:) = []; % remove NaNs
n = size(XY,1);      % number of data points
centroid = mean(XY);   % the centroid of the data set
%     computing moments (note: all moments will be normed, i.e. divided by n)
Mxx=0; Myy=0; Mxy=0; Mxz=0; Myz=0; Mzz=0;
for i=1:n
    Xi = XY(i,1) - centroid(1);  %  centering data
    Yi = XY(i,2) - centroid(2);  %  centering data
    Zi = Xi*Xi + Yi*Yi;
    Mxy = Mxy + Xi*Yi;
    Mxx = Mxx + Xi*Xi;
    Myy = Myy + Yi*Yi;
    Mxz = Mxz + Xi*Zi;
    Myz = Myz + Yi*Zi;
    Mzz = Mzz + Zi*Zi;
end
   
Mxx = Mxx/n;
Myy = Myy/n;
Mxy = Mxy/n;
Mxz = Mxz/n;
Myz = Myz/n;
Mzz = Mzz/n;
%    computing the coefficients of the characteristic polynomial
Mz = Mxx + Myy;
Cov_xy = Mxx*Myy - Mxy*Mxy;
Mxz2 = Mxz*Mxz;
Myz2 = Myz*Myz;
A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz2 - Myz2 - Mz*Mz*Mz;
A0 = Mxz2*Myy + Myz2*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
A22 = A2 + A2;
epsilon=1e-12; 
ynew=1e+20;
IterMax=20;
xnew = 0;
%    Newton's method starting at x=0
for iter=1:IterMax
    yold = ynew;
    ynew = A0 + xnew*(A1 + xnew*(A2 + 4.*xnew*xnew));
    if (abs(ynew)>abs(yold))
        disp('Newton-Pratt goes wrong direction: |ynew| > |yold|');
        xnew = 0;
        break;
    end
    Dy = A1 + xnew*(A22 + 16*xnew*xnew);
    xold = xnew;
    xnew = xold - ynew/Dy;
    if (abs((xnew-xold)/xnew) < epsilon), break, end
    if (iter >= IterMax)
        disp('Newton-Pratt will not converge');
        xnew = 0;
    end
    if (xnew<0.)
        fprintf(1,'Newton-Pratt negative root:  x=%f\n',xnew);
        xnew = 0;
    end
end
%    computing the circle parameters
DET = xnew*xnew - xnew*Mz + Cov_xy;
Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;
Par = [Center+centroid , sqrt(Center*Center'+Mz+2*xnew)];

xCenter = Par(1);
yCenter = Par(2);
radius = Par(3);
end    %    CircleFitByPratt

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
% function [xc,yc,r] = fitcircle(x,y)
% 
% % [xc,yc,r] = fitcircle(x,y)
% %
% % INPUTS:
% %     x, y -- x,y coordinates of a rat running on a circular track
% %
% % OUTPUTS:
% %     xc, yc -- x and y coordinate of the center of the track
% %     r -- radius of the track
% %
% % AT 2011
% 
% xx = x.*x;
% yy = y.*y;
% xy = x.*y;
% 
% A = [nansum(x)  nansum(y)  length(x);
%     nansum(xy) nansum(yy) nansum(y);
%     nansum(xx) nansum(xy) nansum(x)];
% B = [-nansum(xx+yy);
%     -nansum(xx.*y+yy.*y);
%     -nansum(xx.*x+xy.*y)];
% 
% a = A \ B;
% xc = a(1) / -2;
% yc = a(2) / -2;
% r  = sqrt((a(1)^2 + a(2)^2) / 4 - a(3));
% end

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


