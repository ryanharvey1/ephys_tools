function adjustmenu(fitFigure)
%ADJUSTMENU Adjust FitFigure menu
%
%   Helper file for FitFigure will modify some builtin menus and add custom
%   menus.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $    $Date: 2011/05/09 00:40:10 $

sffig = fitFigure.Handle;

h = findall(sffig, 'Type','uimenu', 'Parent',sffig);

% Find the Fit Menu and add more items
h0 = findall(h,'flat', 'Label', 'F&it');
closeFitMenu = uimenu(h0, 'Label', '&Close Fit', ...
    'Tag', 'CloseFitMenuItem', ...
    'Accelerator', 'W', ...
    'Callback', {@iCloseFit, fitFigure});
deleteFitMenu = uimenu(h0, 'Label', '&Delete Fit', ...
    'Tag', 'DeleteFitMenuItem', 'Callback', ...
    {@iDeleteFit, fitFigure});
dupFitMenu = uimenu(h0, 'Label', 'Du&plicate Fit', 'Separator', 'on', ...
    'Tag', 'DuplicateFitMenuItem', 'Callback',  {@iDuplicateFit, fitFigure});

fitFigure.UseDataFromMenu = sftoolgui.DynamicUIMenu(h0, ...
    'Label', '&Use Data From', 'Tag', 'UseDataFromMenuItem');
saveFitToWorkSpaceMenu = uimenu(h0, 'Label', '&Save to Workspace...', ...
    'Separator', 'on', 'Tag', 'SaveToWorkspaceMenuItem', ...
    'Callback', {@iSaveFitToWS, fitFigure});
uimenu(h0, 'Label', 'Specify &Validation Data', ...
    'Tag', 'SpecifyValidationDataMenuItem', ...
    'Separator', 'on', 'Callback',  {@iSpecifyValidationData, fitFigure});

% Create a structure of Fit menu submenus
fitMenus = struct('dupFitMenu', dupFitMenu, ...
    'deleteFitMenu', deleteFitMenu, ...
    'closeFitMenu', closeFitMenu, ...
    'saveFitToWorkSpaceMenu', saveFitToWorkSpaceMenu);

% iSetupFitMenu is the called when the Fit menu is clicked.
set(h0, 'Callback',  {@iSetupFitMenu, fitFigure, fitMenus});

% Add a View Menu
h0 = uimenu(sffig, 'Label', '&View', 'Tag', 'ViewMenu', 'Position', iDesktopMenuPosition(sffig));

fitSettingsMenu = uimenu(h0, 'Label', '&Fit Settings', 'Callback', ...
    {@iToggleFitSettings, fitFigure}, 'Tag', 'SFFitSettingsMenu');
fitResultsMenu = uimenu(h0, 'Label', 'Fit &Results', 'Callback', ...
    {@iToggleFitResults, fitFigure}, 'Tag', 'SFFitResultsMenu');

% Add the plots
plotMenus = cell(length(fitFigure.PlotPanels), 1);
for i=1:length(fitFigure.PlotPanels)
    plotMenus{i} = uimenu(h0, ...
        'Label', fitFigure.PlotPanels{i}.Name,  ...
        'Tag', ['sftool' fitFigure.PlotPanels{i}.Tag 'Menu'], ...
        'Callback', {@iCallViewPlotsMenu, fitFigure, fitFigure.PlotPanels{i}});
end

% set the separator on to the first plot menu item.
set(plotMenus{1}, 'Separator', 'on');

tableOfFitsMenu = uimenu(h0, 'Label', '&Table of Fits', 'Callback', ...
    @iToggleTableOfFits, 'Separator', 'on', 'Tag', 'SFTableOfFitsMenu');

% Create a structure of View menu submenus
viewMenus = struct('fitSettingsMenu', fitSettingsMenu, ...
    'fitResultsMenu', fitResultsMenu, ...
    'tableOfFitsMenu', tableOfFitsMenu);

% iSetupViewMenu is called when the View Menu is clicked
set(h0, 'Callback', {@iSetupViewMenu, fitFigure, viewMenus, plotMenus});

% create the tools menu
toolsMenu = uimenu(sffig, 'Label', '&Tools', 'Tag', 'ToolsMenu', 'Position', iDesktopMenuPosition(sffig));

% if HGUsingMATLABClasses, add zoom in, zoom out, and pan
if curvefit.isHandlegraphics(sffig)
    % Create "empty" menus to pass into the setup method
    zoomInMenu = [];
    zoomOutMenu = [];
    panMenu = [];
else
    zoomInMenu = uimenu(toolsMenu, 'Label', 'Zoom &In', 'Tag', 'exploration.zoom3dInMenu', ...
        'Callback', {@fitFigure.zoomModeCallback, 'in'});
    
    zoomOutMenu = uimenu(toolsMenu, 'Label', 'Zoom &Out', 'Tag', 'exploration.zoom3dOutMenu', ...
        'Callback', {@fitFigure.zoomModeCallback, 'out'});
    
    panMenu = uimenu(toolsMenu, 'Label', '&Pan', 'Tag', 'exploration.pan3dMenu', ...
        'Callback', @fitFigure.panModeCallback);
end

% It is very important that the mode Tag matches the "Data Cursor" tag
% figuretools.m. Menu and toolbar button synchronization depend on this.
dataCursorMenu = uimenu(toolsMenu, 'Label', 'D&ata Cursor', 'Tag', 'figMenuDatatip', ...
    'Callback', 'toolsmenufcn Datatip');

excludeMenu = uimenu(toolsMenu, 'Label', '&Exclude Outliers', 'Tag', ...
    'sftoolExcludeOutlierMenu', 'Callback', @fitFigure.excludeModeCallback);

legendMenu = uimenu(toolsMenu, 'Label', '&Legend', 'Separator', 'on', ...
    'Tag', 'LegendMenu', 'Callback', @fitFigure.toggleLegendState);
gridMenu = uimenu(toolsMenu, 'Label', '&Grid', ...
    'Callback', @fitFigure.toggleGridState, ...
    'Tag', 'GridMenu');

predictionBoundsMenu = uimenu(toolsMenu, 'Label', 'Prediction &Bounds', ...
    'Separator', 'on', 'Tag', 'SurfacePredictionBoundsMenu');
pbmNone = uimenu(predictionBoundsMenu, 'Label', '&None', 'Tag', ...
    'PredBoundsNone', 'Checked', 'on', ...
    'Callback', {@iSetPredictionLevel, fitFigure, 0});
pbm90 = uimenu(predictionBoundsMenu, 'Label', '9&0%', 'Tag', ...
    'PredBounds90', 'Callback', {@iSetPredictionLevel, fitFigure, 90});
pbm95 = uimenu(predictionBoundsMenu, 'Label', '9&5%', 'Tag', ...
    'PredBounds95', 'Callback', {@iSetPredictionLevel, fitFigure, 95});
pbm99 = uimenu(predictionBoundsMenu, 'Label', '9&9%', 'Tag', ...
    'PredBounds99', 'Callback',{@iSetPredictionLevel, fitFigure, 99});
pbmOther = uimenu(predictionBoundsMenu, 'Label', '&Custom...', 'Tag', ...
    'PredBoundsOther', 'Callback', {@iSetCustomPredictionLevel, fitFigure});

% Create a structure of predictionBounds menu submenus
predictionBoundsMenus = struct('pbmNone', pbmNone, ...
    'pbm90', pbm90, ...
    'pbm95', pbm95, ...
    'pbm99', pbm99, ...
    'pbmOther', pbmOther);

% iSetupPredictionBoundsMenu is called when the prediction Bounds menu is
% clicked
set(predictionBoundsMenu, 'Callback', ...
    {@iSetupPredictionBoundsMenu, fitFigure, predictionBoundsMenus});

% Note: unlike other tool menu items, which are enabled/disabled based on
% "no data" selected, the "Axes Limits" item is always enabled.
uimenu(toolsMenu, 'Label', 'A&xes Limits', 'Tag', 'SFAxisLimitControlMenu', ...
    'Separator', 'on', 'Callback', @fitFigure.showAxisLimitsDialog);

% Create a structure of Tools menu submenus
toolsMenus = struct('dataCursorMenu', dataCursorMenu, ...
    'excludeMenu', excludeMenu, ...
    'zoomInMenu', zoomInMenu, ...
    'zoomOutMenu', zoomOutMenu, ...
    'panMenu', panMenu, ...
    'legendMenu', legendMenu, ...
    'gridMenu', gridMenu, ...
    'predictionBoundsMenu', predictionBoundsMenu);

% iSetupToolsMenu is called when the Tools menu is clicked.
set(toolsMenu, 'Callback', {@iSetupToolsMenu, fitFigure, toolsMenus});

% Add a Help Menu
helpMenu = uimenu(sffig, 'Label', '&Help', 'Tag', 'HelpMenu');
uimenu(helpMenu, 'Label', 'Curve &Fitting Tool Help', ...
    'Tag', 'SurfaceFittingToolMenuItem', ...
    'Callback', @(s, e) fitFigure.HSFTool.sftoolHelp());
uimenu(helpMenu, 'Label', '&Curve Fitting Toolbox Help', ...
    'Tag', 'CurveFittingToolHelpMenuItem', ...
    'Callback', @(s, e) fitFigure.HSFTool.cftoolHelp());
uimenu(helpMenu, 'Label', '&Demos', 'Tag', 'DemosHelpMenuItem', ...
    'Callback', @(s, e) fitFigure.HSFTool.demosHelp());
uimenu(helpMenu, 'Label', 'About Curve Fitting &Toolbox', ...
    'Tag', 'AboutCurveFittingMenuItem', ...
    'Callback', @(s, e) fitFigure.HSFTool.aboutHelp());
end

function iSetupToolsMenu(~, ~, fitFigure, toolsMenus)
% function iSetupToolsMenu(src, event, fitFigure, toolsMenus)
% iSetupToolsMenu is called when the tools menu is clicked
dataCursorMenu = toolsMenus.dataCursorMenu;
excludeMenu = toolsMenus.excludeMenu;
zoomInMenu = toolsMenus.zoomInMenu;
zoomOutMenu = toolsMenus.zoomOutMenu;
panMenu = toolsMenus.panMenu;
legendMenu = toolsMenus.legendMenu;
gridMenu = toolsMenus.gridMenu;
predictionBoundsMenu = toolsMenus.predictionBoundsMenu;

if isAnyDataSpecified(fitFigure.HFitdev.FittingData) % Some data is specified
    set( [dataCursorMenu, zoomInMenu, zoomOutMenu, panMenu, predictionBoundsMenu], 'Enable', 'on');
    if isFittingDataValid(fitFigure.HFitdev)
        set(excludeMenu, 'Enable', 'on');
    else
        set(excludeMenu, 'Enable', 'off');
    end
    if fitFigure.LegendOn
        set(legendMenu, 'Checked', 'on', 'Enable', 'on');
    else
        set(legendMenu, 'Checked', 'off', 'Enable', 'on');
    end
    set(gridMenu, 'Checked', fitFigure.GridState, 'Enable', 'on');
    iMakeActiveModeSelected(fitFigure.Handle, toolsMenus);
else % No data specified, disable all the items and make sure checked state is off
    set( [dataCursorMenu, excludeMenu, zoomInMenu, zoomOutMenu, panMenu, legendMenu, gridMenu], 'Enable', 'off', 'Checked', 'off');
    % Prediction bounds doesn't display check, so just disable it.
    set(predictionBoundsMenu, 'Enable', 'off');
end
end

function iMakeActiveModeSelected(hFigure, toolsMenus)
% iMakeActiveModeSelected sets the 'Checked' property to 'on' on the menu
% that corresponds to the active mode. The mutually exclusive modes are:
% rotate3D, data cursor, exclude, zoom in, zoom out and pan. The default
% mode is rotate3D, but there is no menu item or (toolbar button) for that
% mode. It needs to be included here to prevent a warning in the otherwise
% case.

% Get the menus
dataCursorMenu = toolsMenus.dataCursorMenu;
excludeMenu = toolsMenus.excludeMenu;
zoomInMenu = toolsMenus.zoomInMenu;
zoomOutMenu = toolsMenus.zoomOutMenu;
panMenu = toolsMenus.panMenu;

% First make sure that no item is still checked
set( [dataCursorMenu, excludeMenu, zoomInMenu, zoomOutMenu, panMenu], ...
    'Checked', 'off');
hManager = uigetmodemanager(hFigure);
hMode = hManager.CurrentMode;
switch hMode.Name
    case 'Exploration.Rotate3d'
        % No menu item for Rotate3D
    case 'sftoolgui.sfExcludeMode'
        set(excludeMenu, 'Checked', 'on');
    case 'Exploration.Datacursor'
        set(dataCursorMenu, 'Checked', 'on');
    case 'exploration.zoom3d'
        if strcmpi(hMode.ModeStateData.Direction, 'in')
            set(zoomInMenu, 'Checked', 'on');
        else % out
            set(zoomOutMenu, 'Checked', 'on');
        end
    case 'exploration.pan3d'
        set(panMenu, 'Checked', 'on');
    otherwise
        warning(message('curvefit:sftoolgui:FitFigure:adjustmenu:InvalidMode', hMode.Name));
end
end

function iSetupPredictionBoundsMenu(~, ~, fitFigure, pbMenus)
% function iSetupPredictionBoundsMenu(src, event, fitFigure, pbMenus)
% iSetupPredictionBoundsMenu is called when the predictionBoundsMenu is
% clicked.
pbmNone = pbMenus.pbmNone;
pbm90 = pbMenus.pbm90;
pbm95 = pbMenus.pbm95;
pbm99 = pbMenus.pbm99;
pbmOther = pbMenus.pbmOther;

set([pbmNone pbm90 pbm95 pbm99 pbmOther], 'Checked', 'off');
level = fitFigure.HSurfacePanel.PredictionLevel;
switch level
    case 0
        menu = pbmNone;
    case 90
        menu = pbm90;
    case 95
        menu = pbm95;
    case 99
        menu = pbm99;
    otherwise %must be other
        menu = pbmOther;
end
set(menu, 'Checked', 'on');
end

% Tools menu callbacks
function fitFigure = iSetPredictionLevel(~, ~, fitFigure, level)
%iSetPredictionLevel(src, event, fitFigure, level)
fitFigure.HSurfacePanel.PredictionLevel = level;
end

function fitFigure = iSetCustomPredictionLevel(~, ~, fitFigure)
% iSetCustomPredictionBounds(src, event, fitFigure)
levelTxt = inputdlg({xlate('Confidence level (%):')},...
    xlate('Set Confidence Level'), 1, ...
    {fitFigure.OtherPredictionBounds}, 'on');
if ~isempty(levelTxt)
    level = str2double(levelTxt{1});
    if ~isfinite(level) || ~isreal(level) || level<=0 ...
            || level>=100
        errordlg(sprintf(...
            ['Bad confidence level: "%s".\nMust be a percentage ',...
            'larger than 0 and smaller than 100.\n'], levelTxt{1}));
    else
        fitFigure.OtherPredictionBounds = levelTxt{1};
        fitFigure.HSurfacePanel.PredictionLevel =  level;
    end
end
end

function iSetupFitMenu(~, ~, fitFigure, fitMenus)
% function iSetupFitMenu(src, evt, fitFigure, fitMenus)
% iSetupFitMenu is called when the Fit menu is clicked

dupMenu = fitMenus.dupFitMenu;
deleteMenu = fitMenus.deleteFitMenu;
closeMenu = fitMenus.closeFitMenu;
saveFitToWSMenu = fitMenus.saveFitToWorkSpaceMenu;

% Refresh the "Open Fit" menu
refreshOpenFitMenu(fitFigure, @fitFigure.openFit);

% Refresh the "Use Data From" menu
iRefreshUseDataFromMenu(fitFigure, {@iUseDataFrom, fitFigure});

selectedFitName = fitFigure.HFitdev.FitName;

label = sprintf('Du&plicate "%s"', selectedFitName);
set(dupMenu, 'Label', label, 'Enable', 'on');

label = sprintf('&Close "%s"', selectedFitName);
set(closeMenu, 'Label', label, 'Enable', 'on');

label = sprintf('&Delete "%s"', selectedFitName);
set(deleteMenu, 'Label', label, 'Enable', 'on');

% The "Save to Workspace" menu should be enabled if a fitting operation had
% been attempted on the selected fit and there were no errors. Otherwise,
% the menu should be disabled.
if isFitted(fitFigure.HFitdev)
    set(saveFitToWSMenu, 'Enable', 'on');
else
    set(saveFitToWSMenu, 'Enable', 'off');
end
end

%Fit menu callbacks
function iCloseFit(~, ~, fitFigure)
%iCloseFit(src, event, fitFigure)
delete(fitFigure.Handle);
end

function iDeleteFit(~, ~, fitFigure)
%iDeleteFit(src, event, fitFigure)
fitFigure.HSFTool.HFitsManager.deleteFit(fitFigure.FitUUID);
end

function iDuplicateFit(~, ~, fitFigure)
%iDuplicateFit(src, event, fitFigure)
fitFigure.HSFTool.HFitsManager.duplicateFit(fitFigure.FitUUID);
end

function iUseDataFrom(src, ~,  fitFigure)
%iUseDataFrom(src, event, fitFigure)
dataSrcFit = get(src,  'UserData');
fitFigure.HFitdev = updateDataFromADifferentFit(fitFigure.HFitdev, dataSrcFit);
end

function iSaveFitToWS(~, ~, fitFigure)
%iSaveFitToWS(src, event, fitFigure)
fitFigure.HFitdev.saveFitToWorkspace();
end

function iSpecifyValidationData(~, ~, fitFigure)
%iSpecifyValidationData(src, event, fitFigure)
fitFigure.HFittingPanel.showValidationDialog();
end

function iSetupViewMenu(~, ~, fitFigure, viewMenus, plotMenus)
% function iSetupViewMenu(src, event, fitFigure, viewMenus, plotMenus)
% iSetupViewMenu is called when the View menu is clicked. It updates the
% view menu submenus.

% Update the Fit Settings menu
iUpdateFitSettingsMenu(fitFigure, viewMenus.fitSettingsMenu);

% Update the Fit Results menu
iUpdateResultsMenu(fitFigure, viewMenus.fitResultsMenu);

% Update the plot menus
updatePlotControls(fitFigure, plotMenus, 'Checked');

% Update the Table of Fits menu
iUpdateTableOfFitsMenu(viewMenus.tableOfFitsMenu);

end

function iUpdateFitSettingsMenu(fitFigure, fitSettingsMenu)
% iUpdateFitSettingsMenu sets the 'Checked' and 'Enable' properties of the
% Fit Settings menu.
if isAnyDataSpecified(fitFigure.HFitdev.FittingData) % Some data is specified
    % Checked state should match Panel's visible state
    checked = fitFigure.HFittingPanel.Visible;
    enable = 'on';
else % No data specified
    % Set the fit settings menu item to be checked, but disabled.
    checked = 'on';
    enable = 'off';
end
set(fitSettingsMenu, 'Checked', checked, 'Enable', enable);
end

function iUpdateResultsMenu(fitFigure, resultsMenu)
% iUpdateResultsMenu sets the 'Checked' property of the Results menu
set(resultsMenu, 'Checked', fitFigure.HResultsPanel.Visible);
end

function iUpdateTableOfFitsMenu(tableOfFitsMenu)
% iUpdateTableOfFitsMenu sets the 'Checked' property of the Table of Fits
% menu
dt = javaMethodEDT( 'getDesktop', 'com.mathworks.mlservices.MatlabDesktopServices' );
if dt.isClientShowing('Table of Fits')
    checked = 'on';
else
    checked = 'off';
end
set(tableOfFitsMenu, 'Checked', checked);
end

% View menu callbacks
function iToggleFitSettings(src, ~, fitFigure)
%iToggleFitSettings(src, event, fitFigure)
if strcmpi(get(src, 'Checked'), 'on')
    fitFigure.HFittingPanel.Visible = 'off';
else
    fitFigure.HFittingPanel.Visible = 'on';
end
updatePositions(fitFigure);
notify(fitFigure, 'SessionChanged');
end

function iToggleFitResults(src, ~, fitFigure)
%iToggleFitResults(src, event, fitFigure)
if strcmpi(get(src, 'Checked'), 'on')
    fitFigure.HResultsPanel.Visible = 'off';
else
    fitFigure.HResultsPanel.Visible = 'on';
end
updatePositions(fitFigure);
notify(fitFigure, 'SessionChanged');
end

function iCallViewPlotsMenu(src, ~, fitFigure, plotPanel)
% callViewPlots(src, eventdata, fitFigure, plotPanel)
if strcmpi(get(src, 'Checked'), 'on')
    newState = 'off';
else
    newState = 'on';
end
plotPanel.Visible = newState;

if strcmpi(newState, 'on')
    if isa(plotPanel, 'sftoolgui.SurfacePanel')
        plotSurface(fitFigure.HSurfacePanel);
    elseif isa(plotPanel, 'sftoolgui.ResidualsPanel')
        plotResiduals(fitFigure.HResidualsPanel);
    elseif isa(plotPanel, 'sftoolgui.ContourPanel')
        plotSurface(fitFigure.HContourPanel);
    end
end

notify(fitFigure, 'PlotVisibilityStateChanged');
notify(fitFigure, 'SessionChanged');
end

function iToggleTableOfFits(src, ~)
%iToggleTableOfFits(src, event)
dt = javaMethodEDT( 'getDesktop', ...
    'com.mathworks.mlservices.MatlabDesktopServices' );
if strcmpi(get(src, 'Checked'), 'on')
    dt.hideClient('Table of Fits');
else
    dt.showClient('Table of Fits');
end
end

function desktopMenuPosition = iDesktopMenuPosition(fig)
% iDesktopMenuPosition returns the position of the fig's Desktop menu
dm = findall(fig, 'Type', 'uimenu', 'Tag', 'figMenuDesktop'); % Find the desktop menu
desktopMenuPosition = get(dm, 'Position');
end

function iRefreshUseDataFromMenu(fitFigure, callback)
% iRefreshUseDataFromMenu(this, callback) removes existing submenus from
% the "Use Data From" menu and then adds submenus for all existing fits
% (except the currently selected fit) that have valid data.
removeSubmenus(fitFigure.UseDataFromMenu);
HFitdevs = fitFigure.HSFTool.HFitsManager.Fits;
n = size(HFitdevs, 2);
for i = 1:n
    fit = HFitdevs{i};
    % Don't list the selected fit or any fits that don't have valid data
    if ~isAssociatedFit(fitFigure, fit) && isFittingDataValid(fit)
        addSubmenu(fitFigure.UseDataFromMenu, 'Label', fit.FitName, ...
            'Tag', ['UseDataFromFit' fit.FitName 'MenuItem'], ...
            'UserData', fit, 'Callback', callback );
    end
end
end

