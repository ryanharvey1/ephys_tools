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
W0 = sqrt((x - xc).* (x - xc) + (y - yc) .* (y - yc));
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