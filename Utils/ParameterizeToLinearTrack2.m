function [V,W] = ParameterizeToLinearTrack2(X,Y)
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
%
% Updated by Ryan H 6/1/17 (lines 69-72) to stop flipping of data
%
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
% Find placement of nan
%--------------------

ts=1:length(x);
tsnan=ts(isnan(x));
ts_not_nan=ts(~isnan(x));


%--------------------
% Determine eignevectors
%--------------------

Z = [x - nanmean(x)  y - nanmean(y)];
Z=Z(~isnan(Z(:,1)),:);

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

% if (max(V0) < abs(min(V0)))
%   % then we are reversed
%   V0 = -V0;
% end

% zero it
V0 = V0 - min(V0);

V = V0;
W = W0;

% deal with nan values
V=[V,ts_not_nan'];
W=[W,ts_not_nan'];

V=[V;nan(length(tsnan),1),tsnan'];
W=[W;nan(length(tsnan),1),tsnan'];

[~,I]=sort(V(:,2));
V=V(I,1);
W=W(I,1);

end
