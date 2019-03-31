function updateFitTypeValues ( this )
%updateFitTypeValues FitFigure callback to Fitdev's FitTypeFitValuesUpdated event

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $    $Date: 2010/11/01 19:22:12 $

updateResults(this.HResultsPanel, ' ');
clearSurface(this.HSurfacePanel);
clearSurface(this.HContourPanel);
updateInformationPanel(this);
updateLegends(this);
end