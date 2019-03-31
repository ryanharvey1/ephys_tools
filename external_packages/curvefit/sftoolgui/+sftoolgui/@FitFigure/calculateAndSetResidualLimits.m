function calculateAndSetResidualLimits(hFitFigure)
% calculateAndSetResidualLimits calculates and sets the Residual Limits
%
% calculateAndSetResidualLimits (hFitFigure) determines residual limits
% and sets AxesViewModel ResidualLimits property.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.4.1 $    $Date: 2011/07/18 00:31:10 $

setLimit(hFitFigure.AxesViewModel, 'ResidualLimits', iGetResidLimits(hFitFigure));
end

function residLimits = iGetResidLimits(hFitFigure)
% iGetResidLimits gets the limits for residuals
[resids, vResids] = getResiduals(hFitFigure.HResidualsPanel);
if isempty(resids)
    residLimits = [-1 1];
else
    residLimits = cfswitchyard( 'cfAxisLimitsFromData', resids, 0.05);
    if ~isempty(vResids)
        vResidLimits = cfswitchyard( 'cfAxisLimitsFromData', vResids, 0.05);
        residLimits(1) = min(residLimits(1), vResidLimits(1));
        residLimits(2) = max(residLimits(2), vResidLimits(2));
    end
    residLimits = iGetModifiedLimits(hFitFigure, residLimits);
end
end

function limits = iGetModifiedLimits(hFitFigure, limits)
% iGetModifiedLimits adjusts limits to make effective 0 residuals look like
% 0 residual
fittingData = hFitFigure.HFitdev.FittingData;
if isCurveDataSpecified(fittingData) % Curves
    [~, output] = getCurveValues(fittingData);
else % Surface
    [~, ~, output] = getValues(fittingData);
end

if isempty( output )
    range = 1;
else
    range = max( output ) - min( output );
end

minLim = max( 1e-8*range, 1e-5 );

limits(1) = min( limits(1), -minLim );
limits(2) = max( limits(2),  minLim );
end

