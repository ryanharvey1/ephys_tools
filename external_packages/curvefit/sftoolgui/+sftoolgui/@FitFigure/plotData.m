function plotData(hFitFigure, axesNeedPreparation)
%plotData Plot FitFigure fitting and validation data lines
%
% plotData(hFitFigure, axesNeedPreparation) plots hFitFigure fitting and
% validation data lines. axesNeedPreparation specifies whether or not to
% set axes limits and labels based on the data to be plotted. If
% axesNeedPreparation is true, axes limits will be set based on
% the data. If false, axes limits will not be changed.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6.2.1 $    $Date: 2011/07/18 00:31:14 $

hFitdev = hFitFigure.HFitdev;

previewFValues = sftoolgui.util.previewValues(hFitdev.FittingData);
previewVValues = sftoolgui.util.previewValues(hFitdev.ValidationData);

if axesNeedPreparation
    calculateAndSetLimits(hFitFigure);
    prepareAxesForPreview(hFitFigure.HSurfacePanel);
    clearFittingDataLine(hFitFigure.HSurfacePanel);
    clearValidationDataLine(hFitFigure.HSurfacePanel);
    
    prepareAxesForPreview(hFitFigure.HContourPanel);
    clearFittingDataLine(hFitFigure.HContourPanel);
    clearValidationDataLine(hFitFigure.HContourPanel);
    
    prepareAxesForPreview(hFitFigure.HResidualsPanel);
end

% Plot the fitting data if all specified data have the same number of
% elements.
if areNumSpecifiedElementsEqual(hFitdev.FittingData)
    plotDataLineWithExclusions(hFitFigure.HSurfacePanel, ...
        previewFValues, hFitdev.Exclusions, hFitdev.FittingData.Name);
    plotDataLineWithExclusions(hFitFigure.HContourPanel, ...
        previewFValues, hFitdev.Exclusions, hFitdev.FittingData.Name);
    plotDataLineWithExclusions(hFitFigure.HResidualsPanel);
end

% Plot the validation data if all specified data have the same number of
% elements.
if areNumSpecifiedElementsEqual(hFitdev.ValidationData)
    plotDataLine(hFitFigure.HSurfacePanel, ...
        hFitFigure.HSurfacePanel.HValidationDataLine, previewVValues, hFitdev.ValidationData.Name);
    plotDataLine(hFitFigure.HContourPanel, ...
        hFitFigure.HContourPanel.HValidationDataLine, previewVValues, hFitdev.ValidationData.Name);
end

updateLegends(hFitFigure)
end
