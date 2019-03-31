function createToolbar(fitFigure)
%createToolbar Create FitFigure toolbar
%
%   createToolbar - helper file to create a toolbar for FitFigures

%   Copyright 2009-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $    $Date: 2011/05/09 00:40:11 $

% Create the toolbar
hToolbar = uitoolbar('Parent', fitFigure.Handle);

% Create plot panel buttons
iCreatePlotPanelButtons(fitFigure, hToolbar);

% if HGUsingMATLABClasses, add zoom in, zoom out, and pan
if ~curvefit.isHandlegraphics(fitFigure.Handle)
    % Create zoom in button
    iCreateZoomInButton(fitFigure, hToolbar);
    
    % Create zoom out button
    iCreateZoomOutButton(fitFigure, hToolbar);
    
    % Create pan button
    iCreatePanButton(fitFigure, hToolbar);
end

% Create data cursor button
iCreateDataCursorButton(fitFigure, hToolbar);

% Create exclude button
iCreateExcludeButton(fitFigure, hToolbar);

% Create legend button
iCreateLegendButton(fitFigure, hToolbar);

% Create grid button
iCreateGridButton(fitFigure, hToolbar)

% Create axis limits dialog button
iCreateAxisLimitsDialogButton(fitFigure, hToolbar)

end

% Data Cursor toolbar button
function iCreateDataCursorButton(fitFigure, hToolbar)
dataCursorButton = uitoolfactory(hToolbar, 'Exploration.DataCursor');
set(dataCursorButton, 'TooltipString', xlate('Data cursor'));
iSetDataListener(fitFigure, dataCursorButton);
iEnableButton(fitFigure, dataCursorButton);
% if ~HGUsingMATLABClasses, add a separator
if curvefit.isHandlegraphics(fitFigure.Handle)
    set(dataCursorButton, 'Separator', 'on');
end
end

% Zoom in toolbar button
function iCreateZoomInButton(fitFigure, hToolbar)
% Tag needs to match sfzoom3d designation
zoomInButton = uitoggletool(hToolbar, ...
    'TooltipString', xlate('Zoom in'), ...
    'ClickedCallback', {@fitFigure.zoomModeCallback, 'in'}, ...
    'Tag', 'exploration.zoom3din');
% Add a separator
set(zoomInButton, 'Separator', 'on');
iSetMatlabIcon(zoomInButton, 'tool_zoom_in.png');
iSetDataListener(fitFigure, zoomInButton);
iEnableButton(fitFigure, zoomInButton);
end

% Zoom out toolbar button
function iCreateZoomOutButton(fitFigure, hToolbar)
% Tag needs to match sfzoom3d designation
zoomOutButton = uitoggletool(hToolbar, ...
    'TooltipString', xlate('Zoom out'), ...
    'ClickedCallback', {@fitFigure.zoomModeCallback, 'out'}, ...
    'Tag', 'exploration.zoom3dout');
iSetMatlabIcon(zoomOutButton, 'tool_zoom_out.png');
iSetDataListener(fitFigure, zoomOutButton);
iEnableButton(fitFigure, zoomOutButton);
end

% Pan toolbar button
function iCreatePanButton(fitFigure, hToolbar)
% Tag needs to match sfpan3d designation
panButton = uitoggletool(hToolbar, ...
    'TooltipString', xlate('Pan'), ...
    'ClickedCallback', @fitFigure.panModeCallback, ...
    'Tag', 'exploration.pan3d');
iSetMatlabIcon(panButton, 'tool_hand.png');
iSetDataListener(fitFigure, panButton);
iEnableButton(fitFigure, panButton);
end

% Exclude toolbar button
function iCreateExcludeButton(fitFigure, hToolbar)
excludeButton = uitoggletool(hToolbar, ...
    'TooltipString', xlate('Exclude outliers'), ...
    'ClickedCallback', @fitFigure.excludeModeCallback, ...
    'Tag', 'sftoolExcludeOutliersToolbarButton');
sftoolgui.setIcon(excludeButton, sftoolgui.iconPath('excludeSMALL.png'));
% setDataListener
set(excludeButton, 'UserData', event.listener(fitFigure.HFitdev, ...
    'FittingDataUpdated', @(s, e)iUpdateExcludeToolbar(fitFigure, excludeButton)));
% Make sure button state is properly initialized
iUpdateExcludeToolbar(fitFigure, excludeButton);
end

function iUpdateExcludeToolbar(fitFigure, excludeButton)
% iUpdateExcludeToolbar enables/disables the exclude button depending of
% the validity of the Fitting data. It differs from other tool and view
% buttons that merely check for any data specified.
if isFittingDataValid(fitFigure.HFitdev)
    set(excludeButton, 'Enable', 'on');
else
    set(excludeButton, 'Enable', 'off');
end
end

% Plot panels toolbar buttons
function iCreatePlotPanelButtons(fitFigure, hToolbar)
% Set up array to store plot buttons
plotButtons = cell(length(fitFigure.PlotPanels), 1);
% Create a button for each plot panel
for i=1:length(fitFigure.PlotPanels)
    plotPanel = fitFigure.PlotPanels{i};
    plotButtons{i} = uitoggletool(hToolbar, ...
        'State', plotPanel.Visible, ...
        'TooltipString', plotPanel.Description,  ...
        'ClickedCallback', ...
        {@iTogglePlotPanelVisibility, fitFigure, plotPanel}, ...
        'Tag', ['sftool' plotPanel.Tag 'ToolbarButton']);
    sftoolgui.setIcon(plotButtons{i}, sftoolgui.iconPath(plotPanel.Icon));
end

% Create a listeners for PlotVisibilityStateChanged and FittingDataUpdated events
listeners = {event.listener(fitFigure, ...
    'PlotVisibilityStateChanged', @(s, e) updatePlotControls(fitFigure, plotButtons, 'State')), ...
    event.listener(fitFigure.HFitdev, ...
    'FittingDataUpdated', @(s, e) updatePlotControls(fitFigure, plotButtons, 'State')) };

% Store listeners in the first button's user data
set(plotButtons{1}, 'UserData', listeners);

% Make sure button states are properly initialized
updatePlotControls(fitFigure, plotButtons, 'State');
end

function iTogglePlotPanelVisibility(src, ~, fitFigure, plotPanel)
%iTogglePlotPanelVisibility(src, event, fitFigure, plotPanel)
% togglePotPanelVisibility is called when a plot panel toolbar button is
% clicked.
newState = get(src, 'State');
plotPanel.Visible = newState;
if strcmpi(newState, 'on')
    if isa(plotPanel, 'sftoolgui.SurfacePanel')
        plotSurface(plotPanel);
    elseif isa(plotPanel, 'sftoolgui.ResidualsPanel')
        plotResiduals(plotPanel);
    elseif isa(plotPanel, 'sftoolgui.ContourPanel')
        plotSurface(plotPanel)
    end
end
notify(fitFigure, 'PlotVisibilityStateChanged');
notify(fitFigure, 'SessionChanged');
end

% Legend toolbar button
function iCreateLegendButton(fitFigure, hToolbar)
% Create our own legend button rather than using:
% "legend = uitoolfactory(hToolbar, 'Annotation.InsertLegend')"
% because we want the button to match the menu item which reflects (and
% sets) the "LegendOn" property.

legendButton = uitoggletool(hToolbar, ...
    'TooltipString', xlate('Legend'),  ...
    'ClickedCallback', @fitFigure.toggleLegendState, ...
    'Tag', 'sftoolLegendToolbarButton', ...
    'Separator', 'on');
sftoolgui.setIcon(legendButton, sftoolgui.iconPath('tool_legend.png', 'matlab'));

% Create a listener for the LegendStateChanged and FittingDataUpdated
% events
listeners = {event.listener(fitFigure, 'LegendStateChanged', ...
    @(s, e) iUpdateLegendToolbar(fitFigure, legendButton) ), ...
    event.listener(fitFigure.HFitdev, 'FittingDataUpdated', ...
    @(s, e) iUpdateLegendToolbar(fitFigure, legendButton) )};

% Store the listeners in button's UserData
set(legendButton, 'UserData', listeners)

% Make sure button states are properly initialized
iUpdateLegendToolbar(fitFigure, legendButton);
end

function iUpdateLegendToolbar(fitFigure, legendButton)
% Event listener Callback used to update the legend toolbar
% button's State property when the LegendStateChanged event is
% fired.
if isAnyDataSpecified(fitFigure.HFitdev.FittingData)
    if fitFigure.LegendOn
        set(legendButton, 'State', 'on', 'Enable', 'on')
    else
        set(legendButton, 'State', 'off', 'Enable', 'on')
    end
else
    set(legendButton, 'State', 'off', 'Enable', 'off');
end
end

% Grid toolbar button
function iCreateGridButton(fitFigure, hToolbar)
% Toolbar button: grid
gridButton = uitoggletool(hToolbar, ...
    'State', fitFigure.GridState,...
    'TooltipString', xlate('Grid'),...
    'ClickedCallback', @fitFigure.toggleGridState,...
    'Tag','sftoolGridToolbarButton');
iSetCftoolIcon(gridButton, 'grid');

% Create a listener for the GridStateChanged and FittingDataUpdated events
listeners = {event.listener(fitFigure, 'GridStateChanged', ...
    @(s, e) iUpdateGridToolbar(fitFigure, gridButton) ), ...
    event.listener(fitFigure.HFitdev, ...
    'FittingDataUpdated', @(s, e)iUpdateGridToolbar(fitFigure, gridButton))};

% Store the listeners in the button's user data
set(gridButton, 'UserData', listeners);

% Make sure button states are properly initialized
iUpdateGridToolbar(fitFigure, gridButton)
end

function iUpdateGridToolbar(fitFigure, gridButton)
% Event listener Callback used to update the grid toolbar button's
% State property when the GridStateChanged event is fired or when data
% changes.
if isAnyDataSpecified(fitFigure.HFitdev.FittingData)
    set(gridButton, 'Enable', 'on', 'State', fitFigure.GridState);
else
    set(gridButton, 'Enable', 'off', 'State', 'off');
end
end

% Axis Limits Dialog toolbar button
function iCreateAxisLimitsDialogButton(fitFigure, hToolbar)
% Toolbar button: Axis limits dialog
% Note: unlike many of the other toolbar buttons, the axis limits button is
% NOT enabled/disabled depending on "no data selected". It should always be
% enabled.
axisButton = uipushtool(hToolbar, ...
    'TooltipString', xlate('Adjust axes limits'), ...
    'Separator','on',...
    'ClickedCallback',@fitFigure.showAxisLimitsDialog, ...
    'Tag', 'sftoolAxisLimitsDialogToolbarButton');
sftoolgui.setIcon(axisButton, sftoolgui.iconPath('axisLimitDlg.png'));
end

function iSetDataListener(fitFigure, button)
% iSetDataListener creates a listener to the 'FittingDataUpdated' event and
% stores it in the buttons UserData
set(button, 'UserData', event.listener(fitFigure.HFitdev, ...
    'FittingDataUpdated', @(s, e)iEnableButton(fitFigure, button)));
end

function iEnableButton(fitFigure, button)
% iEnableButton enables/disables a button depending on whether or not data
% has been specified
if isAnyDataSpecified(fitFigure.HFitdev.FittingData)
    set(button, 'Enable', 'on');
else
    set(button, 'Enable', 'off');
end
end

function iSetMatlabIcon(button, filename)
% iSetMatlabIcon sets a button's CData with values found in a MATLAB icons
% file
ICONROOT = fullfile(matlabroot,'toolbox','matlab','icons',filesep);
[cdata, ~, alpha] = imread([ICONROOT, filename],'Background','none');
% Converting 16-bit integer colors to MATLAB colorspec
cdata = double(cdata) / 65535.0;
% Set all transparent pixels to be transparent (nan)
cdata(alpha==0) = NaN;
set(button, 'CData', cdata);
end

function iSetCftoolIcon(button, icon)
% sets a button's CData with value found in the cficons.mat file
cficonsFile = fullfile(matlabroot, 'toolbox', 'curvefit', ...
    'cftoolgui', 'private', 'cficons.mat');
if exist(cficonsFile,'file')~=2
    error(message('curvefit:FitFigure:MissingIconFile', 'cficons.mat'));
end
icons = load(cficonsFile,'icons');
set(button, 'CData', icons.icons.(icon));
end
