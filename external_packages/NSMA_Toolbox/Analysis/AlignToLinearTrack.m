function [V,W] = AlignToLinearTrack(X,Y,x0,y0,x1,y1)
%
% [V,W] = AlignToLinearTrack(X,Y,x0,y0,x1,y1)
%
% INPUTS:
%    X,Y -- ctsd or tsd or straight arrays of data
%    x0,y0 & x1,y1 - two points along track (usually ends of the line)
%
% OUTPUTS:
%    V, W -- ctsd or tsd or straight arrays of data
%    V = position along track, W = orthogonal position
%    output is of same type as input

% ADR 1998
%  version L4.0
%  status: PROMOTED

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

vx = x1 - x0; 
vy = y1 - y0;
lenV = sqrt(vx * vx + vy * vy);
vx = vx/lenV; 
vy = vy/lenV;

if y1 == x0
  wx = y1 + 10 - x0; 
  wy = x1 - y0;
else
  wx = y1 - x0; 
  wy = x1 - y0;
end

lenT = wx * vx + wy * vy;
wx = wx - lenT/lenV * (x1 - x0);
wy = wy - lenT/lenV * (y1 - y0);
lenW = sqrt(wx * wx + wy * wy);
wx = wx/lenW; 
wy = wy/lenW;

%--------------------
% Project
%--------------------

V0 = (x - x0) * vx + (y - y0) * vy;
W0 = (x - x0) * wx + (y - y0) * wy;

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