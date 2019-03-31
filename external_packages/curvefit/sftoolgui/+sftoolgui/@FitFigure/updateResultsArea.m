function updateResultsArea(this)
%updateResultsArea FitFigure utility
%
%   updateResultsArea will update the FitFigure's Results area

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/08/29 08:22:35 $

fitObject = this.HFitdev.Fit;
updateInformationPanel(this);

% fitObject will be empty if we're duplicating an "incomplete" fit
if ~isempty(fitObject)
    updateResults(this.HResultsPanel, char( this.HFitdev.ResultsTxt ) );
end
end
