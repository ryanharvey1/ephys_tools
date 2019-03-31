function [V,W] = ParameterizeToCircularTrack(X, Y)

% [V, W] = ParameterizeToCircularTrack(X,Y)
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

%------------------------
% first calculate xc, yc, x0, y0
% use means for center, first point for 0

switch class(X)
   case {'tsd','ctsd'}
      xc = mean(Data(X));
      x0 = Data(X, StartTime(X));
   otherwise 
      xc = mean(X);
      x0 = X(1);
end
   
switch class(Y)
   case {'tsd','ctsd'}
      yc = mean(Data(Y));
      y0 = Data(Y, StartTime(Y));
   otherwise 
      yc = mean(Y);
      y0 = Y(1);
end

[V,W] = AlignToCircularTrack(X,Y, xc, yc, x0, y0);
