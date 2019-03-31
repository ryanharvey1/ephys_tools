function plotFitAndResids(this)
%plotFitAndResids Update FitFigure plot panels and their positions

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.2.2.1 $    $Date: 2011/07/18 00:31:15 $

% calculate and set limits based on fit
calculateAndSetResidualLimits(this);
if strcmp(this.HSurfacePanel.Visible, 'on')
    plotSurface(this.HSurfacePanel);
end
if strcmp(this.HContourPanel.Visible, 'on')
    plotSurface(this.HContourPanel);
end
if strcmp(this.HResidualsPanel.Visible, 'on')
    plotResiduals(this.HResidualsPanel);
end
end