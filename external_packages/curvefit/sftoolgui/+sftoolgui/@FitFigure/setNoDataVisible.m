function setNoDataVisible( this )
%setNoDataVisible sets the "No Data" and Plot Panels Visible state.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $    $Date: 2011/04/11 16:09:40 $

% If data is specified, the plot panels should be visible and the
% "No Data" panel should be invisible.
if isAnyDataSpecified(this.HFitdev.FittingData)
    this.HNoDataPanel.Visible = 'off';
    this.HPlotPanel.Visible = 'on';
else
    % Otherwise, the plot panels should be visible and the "No Data"
    % panel should be invisible.
    this.HPlotPanel.Visible = 'off';
    this.HNoDataPanel.Visible = 'on';
end
end
