function [V,W] = ManualParamToCircularTrack(X, Y)

% [V, W] = ManualParamToCircularTrack(X,Y)
%
% INPUTS:
%     X, Y -- ctsd or tsd or straight arrays of data
%
% OUTPUTS:
%     V, W -- ctsd or tsd or straight arrays of data
%     V = position along circle, W = radial distance from center
%     output is of same type as input

% ADR 1998
%  version L4.0
%  status: PROMOTED
% ZN 2010: asks user to click on starting position on track (or reversal 
%   point of horseshoe track) and uses that as x0,y0
%   (ParameterizeToCircularTrack used the first location for this value
% calls AlignToCircularTrack

%------------------------
% first calculate xc, yc, x0, y0
% asks user to click on center of maze and starting point

figure; plot(Data(X), Data(Y), 'k-')
title('Click on center of track', 'color', 'r')
xlabel('X')
ylabel('Y')
axis equal
disp('Click center of circular track')
[xc,yc] = ginput(1);
title('Click on starting point on track (turn-around point of horseshoe track)', 'color', 'b')
disp('Click on rats starting point on track (turn around point of horseshoe track)')
[x0,y0] = ginput(1);

% switch class(X)
%    case {'tsd','ctsd'}
%       xc = (max(Data(X))+min(Data(X)))/2;
%       %x0 = Data(X, StartTime(X));
%    otherwise 
%       xc = (max(X)+min(X))/2;
%       %x0 = X(1);
% end
%    
% switch class(Y)
%    case {'tsd','ctsd'}
%       yc = (max(Data(Y))+min(Data(Y)))/2;
%       %y0 = Data(Y, StartTime(Y));
%    otherwise 
%       yc = (max(Y)+min(Y))/2;
%       %y0 = Y(1);
% end

[V,W] = AlignToCircularTrack(X,Y, xc, yc, x0, y0);
