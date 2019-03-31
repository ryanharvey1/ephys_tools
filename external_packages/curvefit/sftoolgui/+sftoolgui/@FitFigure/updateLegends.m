function updateLegends(this)
%updateLegends Update all FitFigure plot panels' legends
%
%   updateLegends is called to update all FitFigure plot panels' legends

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/08/29 08:22:33 $

updateLegend(this.HResidualsPanel, this.LegendOn);
updateLegend(this.HSurfacePanel, this.LegendOn);
updateLegend(this.HContourPanel, this.LegendOn);
end