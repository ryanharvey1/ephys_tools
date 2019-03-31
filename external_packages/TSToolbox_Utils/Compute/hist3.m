function [h X Y] = hist3(Data, XBins, YBins, ZBins)
% out = hist2(Data, XBins, YBins)
%
% Makes a 3d histogram of the data.
% XBins and YBins are optional arguments
% which give the number of grid segments.
% - default is 50.

if length(XBins) == 1
    MinX = min(Data(:,1));
    MaxX = max(Data(:,1));
else
    MinX = XBins(1);
    MaxX = XBins(end);
    XBins = length(XBins);
end
if length(YBins) == 1
   MinY = min(Data(:,2));
   MaxY = max(Data(:,2));
else
    MinY = YBins(1);
    MaxY = YBins(end);
    YBins = length(YBins);
end
if length(ZBins) == 1
   MinZ = min(Data(:,3));
   MaxZ = max(Data(:,3));
else
    MinZ = ZBins(1);
    MaxZ = ZBins(end);
    ZBins = length(ZBins);
end

XBin = floor(1 + XBins*(Data(:,1) - MinX) / (MaxX - MinX));
YBin = floor(1 + YBins*(Data(:,2) - MinY) / (MaxY - MinY));
ZBin = floor(1 + ZBins*(Data(:,3) - MinZ) / (MaxZ - MinZ));

ix = isnan(XBin) | isnan(YBin) |  isnan(ZBin);
XBin(ix) = [];
YBin(ix) = [];
ZBin(ix) = [];

h = zeros(XBins, YBins, ZBins);

XBin(find(XBin == XBins+1)) = XBins;
YBin(find(YBin == YBins+1)) = YBins;
ZBin(find(ZBin == ZBins+1)) = ZBins;

for i = 1:size(XBin,1)
    try
    h(XBin(i), YBin(i), ZBin(i)) = h(XBin(i), YBin(i), ZBin(i)) + 1;
    catch
        keyboard
    end
end;

X = [MinX:(MaxX-MinX)/(XBins-1):1.000000001*MaxX];
Y = [MinY:(MaxY-MinY)/(YBins-1):1.000000001*MaxY];
Z = [MinZ:(MaxZ-MinZ)/(ZBins-1):1.000000001*MaxZ];


%  imagesc(h)

% Written by Kenneth D. Harris, adapted by Adrien Peyrache
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu