function calculateAndSetLimits(hFitFigure)
% calculateAndSetLimits determines and sets AxesViewModel limit properties
%
% calculateAndSetLimits(hFitFigure) determines limits based on data and
% residuals and sets the AxesViewModel limit properties.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.4.1 $    $Date: 2011/07/18 00:31:09 $

hFitdev = hFitFigure.HFitdev;

% Get preview values for fitting data
fValues = sftoolgui.util.previewValues(hFitdev.FittingData);
% Get preview values for validation data
vValues = sftoolgui.util.previewValues(hFitdev.ValidationData);

% Get the limits
xlim = iGetLimits(fValues{1}, vValues{1});
ylim = iGetLimits(fValues{2}, vValues{2});
zlim = iGetLimits(fValues{3}, vValues{3});

% Get residual limits.
residualLim = iGetResidLimits(hFitFigure);
if isCurveDataSpecified( hFitFigure.HFitdev.FittingData )
    setLimits(hFitFigure.AxesViewModel, {xlim}, ylim, residualLim);
else
    setLimits(hFitFigure.AxesViewModel, {xlim, ylim}, zlim, residualLim);
end
end

function limits = iGetLimits(fValue, vValue)
% iGetLimits gets limits using the "plotting" rather than real data.
% "plotting" data contains zeros for undefined data values.
if isempty(fValue) && isempty(vValue)
    limits = [-1 1];
elseif isempty(fValue) % only validation value is available
    limits = cfswitchyard( 'cfAxisLimitsFromData', vValue, 0.05 );
elseif isempty(vValue) % only fitting value is specified
    limits = cfswitchyard( 'cfAxisLimitsFromData', fValue,  0.05  );
else % both fitting value and validation value are available;
    % adjust limits to account for both
    limitsFV = cfswitchyard('cfAxisLimitsFromData', fValue,  0.05 );
    limitsVV = cfswitchyard('cfAxisLimitsFromData', vValue,  0.05 );
    limits(1) = min(limitsFV(1), limitsVV(1));
    limits(2) = max(limitsFV(2), limitsVV(2));
end
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

