function updateValidationData( this )
%updateValidationData FitFigure callback to Fitdev's ValidationDataUpdated event 

%   Copyright 2008-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $    $Date: 2010/03/01 05:19:49 $

hFitdev = this.HFitdev;
updateValidationMessage(this.HFittingPanel, getMessageString(hFitdev.ValidationData), getMessageLevel(hFitdev.ValidationData));
plotData(this, true);
updateResultsArea(this);
plotValidationData(this.HResidualsPanel);
end
