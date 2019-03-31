function [V,W] = ParameterizeToLinearTrack(X,Y)
%
% [V,W] = ParameterizeToLinearTrack(X,Y)
% 
% INPUTS:
%    X,Y -- ctsd or tsd or straight arrays of data
%
% OUTPUTS:
%    V,W -- ctsd or tsd or straight arrays of data
%    V = position along track, W = orthogonal position
%    output is of same type as input
% 
% ALGO:
%    Uses SVD to determine linearization.

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
% Determine eignevectors
%--------------------

Z = [x - mean(x)  y - mean(y)];
cv = [Z' * Z];
[v,d] = eig(cv);

A = Z * v(:,1);
B = Z * v(:,2);

%--------------------
% Which is the long trajectory?
%--------------------

if (d(1,1) > d(2,2))
  V0 = A;
  W0 = B;
  veig=v(:,1);
  weig=v(:,2);
else
  V0 = B;
  W0 = A;
  veig=v(:,2);
  weig=v(:,1);
end

if (max(V0) < abs(min(V0)))
  % then we are reversed
  V0 = -V0;
end

% zero it
V0 = V0 - min(V0);

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