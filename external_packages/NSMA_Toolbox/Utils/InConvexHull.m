function b = InConvexHull(px, py, cx, cy, varargin)

% InConvexHull  Determines whether specified points are within a specified convex hull
%
% b = InConvexHull(px, py, cx, cy, varargin)
%
% INPUTS:
%       px, py - list of points to be checked
%       cx, cy - x and y coordinates for a convex hull
%       varargin PARAMETERS:
%           debug (1/0) - default 0
% OUPUTS:
%       b - vector of same length as px & py - 
%           contains 1's where corresponding points are in convex hull, 0's otherwise
%
% ALGORITHM
%   for each point, select a point Q that is definitely outside the polygon
%   in this case Qx = Px & Qy > max(Py)
%   count the number of lines that cross the PQ line
%   if the number of lines is odd then point P is in the polygon
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status PROMOTED


debug = 0;
Extract_varargin;

nConvHullPoints = length(cx);
if length(cy) ~= nConvHullPoints
   error('length(cx) does not match length(py)');
end
nPoints = length(px);
if length(py) ~= nPoints
   error('lengeth(px) does not match length(py)');
end

% start by only including points in xrange
possPoints = find(px > min(cx) & px < max(cx));

if (debug)
   newFig = figure;
end

b = zeros(nPoints,1);
for iK = 1:(nConvHullPoints-1)
   % Does the line PQ intersect the line AB?
   
   Ax = cx(iK); Ay = cy(iK);
   Bx = cx(iK+1); By = cy(iK+1);
   if (Ax > Bx)
      tmp = Ax; Ax = Bx; Bx = tmp;
      tmp = Ay; Ay = By; By = tmp;
      clear tmp
   end
   
   if (Ax ~= Bx)
      intersecty = Ay + (px(possPoints) - Ax)/(Bx - Ax) * (By - Ay);  
      b(possPoints) = b(possPoints) + (Ax < px(possPoints) & Bx > px(possPoints) & py(possPoints) > intersecty);
   end
   
   if (debug)
      figure(newFig); clf; hold on
      f = find(b == 0);
      plot(px(f), py(f), 'b.');
      f = find(b == 1);
      plot(px(f), py(f), 'g.');
      f = find(b == 2);
      plot(px(f), py(f), 'r.');
      plot(cx([iK iK+1]), cy([iK iK+1]), 'm-');
   end
end

if (debug)
   close(newFig)
end
b = rem(b,2);