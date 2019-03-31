function updateFitFigure(this, configuration)
%updateFitFigure Update FitFigure with configuration information
%
%   updateFitFigure(THIS, CONFIGURATION) will update the given FitFigure,
%   THIS, with the given CONFIGURATION information.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $    $Date: 2011/05/09 00:40:19 $

this.LegendOn = configuration.LegendOn;
this.GridState = configuration.Grid;

this.HResidualsPanel.Visible = configuration.ResidualsConfig.Visible;
this.HContourPanel.Visible = configuration.ContourConfig.Visible;
this.HSurfacePanel.Visible = configuration.SurfaceConfig.Visible;
this.HSurfacePanel.PredictionLevel = configuration.SurfaceConfig.PredictionLevel;

plotFitAndResids(this);

iUpdateFittingPanel(this, configuration);

this.HResultsPanel.Visible = configuration.ResultsConfig.Visible;

updateLegends(this);
updateGrids(this);
updatePositions(this);
% send notifications so toolbar states can be set properly
notify(this, 'PlotVisibilityStateChanged');
end

function iUpdateFittingPanel(this, configuration)
% iUpdateFittingPanel updates the fitting panel's 'Visible' property.

% If no data is specified, override the configuration information to make
% sure the fitting panel is visible.
if ~isAnyDataSpecified(this.HFitdev.FittingData)
    this.HFittingPanel.Visible = 'on';
else
    this.HFittingPanel.Visible = configuration.FittingConfig.Visible;
end
end