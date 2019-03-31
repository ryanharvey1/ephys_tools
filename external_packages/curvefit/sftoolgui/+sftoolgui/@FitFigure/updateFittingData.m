function updateFittingData( this, axesNeedPreparation)
%updateFittingData updates plots and results when fitting data changes.
%
% updateFittingData also turns exclude mode off if data is not valid.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $    $Date: 2011/04/11 16:09:42 $

updateResults(this.HResultsPanel, ' ');
updateInformationPanel(this);
clearSurface(this.HSurfacePanel);
clearSurface(this.HContourPanel);
plotData(this, axesNeedPreparation);
% If fitting data is not valid, make sure exclusion mode is
% "off".
if ~isFittingDataValid(this.HFitdev)
    sftoolgui.sfExcludeMode(this, 'off');
end
end
