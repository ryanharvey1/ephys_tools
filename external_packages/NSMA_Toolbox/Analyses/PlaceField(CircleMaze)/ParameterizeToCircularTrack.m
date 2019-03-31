function [V,W,xc,yc] = ParameterizeToCircularTrack(X,Y)

% ParameterizeToCircularTrack  Reparameterizes x,y position data to degrees on circle
%
% [V,W,xc,yc] = ParameterizeToCircle_Shigenori(X,Y)
%
% INPUTS:
%   X,Y - (c)tsd or straight arrays of x,y coordinates
% OUTPUTS:
%   V,W - (c)tsd or straight arrays of data, depending on form of inputs:
%       V = position along circle (degrees), W = radial distance from center
%   xc,yc = center x,y coordinates
%
% First calculate xc, yc, x0, y0
% use means for center, first point for 0
%
% ADR 1998, last modified '04 by MN


switch class(X)
    case {'tsd','ctsd'}
        XData = Data(X);
    otherwise
        XData = X;    
end
       
switch class(Y)
    case {'tsd','ctsd'}
        YData = Data(Y);
    otherwise
        YData = Y;    
end

LRX = XData(find(abs(YData-mean(YData))<10));
LX = LRX(find(LRX < mean(XData)));
RX = LRX(find(LRX > mean(XData)));
UDY = YData(find(abs(XData-mean(XData))<10));
UY = UDY(find(UDY > mean(YData)));
DY = UDY(find(UDY < mean(YData)));

xc = mean([mean(LX) mean(RX)]); % center coordinates
yc = mean([mean(UY) mean(DY)]);

x0 = XData(1);  % starting coordinates
y0 = YData(1);
   

[V,W] = AlignToCircularTrack(X,Y, xc, yc);
