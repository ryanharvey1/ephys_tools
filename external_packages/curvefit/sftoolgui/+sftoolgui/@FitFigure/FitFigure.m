classdef FitFigure < sftoolgui.Figure
    %FitFigure Surface Fitting Tool Fit Figure
    
    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.10.2.1 $    $Date: 2011/07/18 00:31:08 $
    
    events
        SessionChanged
    end
    
    events (NotifyAccess = private)
        %GridStateChanged-- fired when GridState property is changed
        GridStateChanged
        %LegendStateChangedEvent-- fired when LegendOn property is changed
        LegendStateChanged
        %PlotVisibilityStateChanged-- fired when plot panels Visible
        %properties are changed
        PlotVisibilityStateChanged
    end
    
    properties
        HFitdev ;
        FitUUID ;
        HSurfacePanel ;
        HResidualsPanel ;
        HFittingPanel ;
        HResultsPanel ;
        HContourPanel ;
        OtherPredictionBounds = '95';
        HLimitSpinnersDialog ;
        AxesViewModel ;
    end
    
    properties(AbortSet)
        LegendOn = true;
        GridState = 'on';
    end
    
    properties (Dependent = true, SetAccess = private)
        Configuration ;
    end
    
    properties(SetAccess = 'private')
        HNoDataPanel ;
        HPlotPanel ;
    end
    
    properties(SetAccess = 'private', GetAccess = 'private')
        FitdevListeners = {};
        PlotPanels ;
        UseDataFromMenu ;
    end
    
    methods
        function this = FitFigure(sftool, hFitdev, config)
            %FITFIGURE Constructor for the surface fitting tool fit
            %figure
            %
            %   H = FITFIGURE
            this = this@sftoolgui.Figure(sftool);
            
            % create an AxesViewModel
            this.AxesViewModel = sftoolgui.AxesViewModel();
            
            %set the Figure color to be the same as the panels
            uipanelColor = get(this.Handle, 'DefaultUipanelBackgroundColor');
            set(this.Handle, 'Color', uipanelColor);
            
            this.HFitdev = hFitdev;
            this.FitUUID = hFitdev.FitID;
            
            % Create listeners to Fitdev
            this.FitdevListeners = {
                hFitdev.addlistener( 'FitTypeFitValuesUpdated', @(s, e) this.updateFitTypeValues()  )
                hFitdev.addlistener( 'FittingDataUpdated',      @(s, e) this.updateFittingData(true)    )
                hFitdev.addlistener( 'FittingDataUpdated',      @(s, e) this.setNoDataVisible()    )
                hFitdev.addlistener( 'ValidationDataUpdated',   @(s, e) this.updateValidationData() )
                hFitdev.addlistener( 'FitFitted',               @(s, e) this.fitFitted()            )
                hFitdev.addlistener( 'FitNameUpdated',          @(s, e) this.fitNameUpdated()       )
                hFitdev.addlistener( 'ExclusionsUpdated',       @(s, e) this.updateFittingData(false)   )
                };
            
            % Note that for the ExclusionsUpdated event, i.e., when exclusions
            % have been updated, we want to do the same things as when fitting
            % data has been updated
            
            % Set DeleteFcn
            set(this.Handle, 'DeleteFcn', @this.fitFigureDeleteFcn);
            
            % Create panels
            this = iCreatePanels(this);
            
            iMakeFigureGood(this);
            setNoDataVisible(this);
            
            updateFittingData(this, true);
            if ~isempty(config)
                updateFitFigure(this, config);
            end
            if isFitted(this.HFitdev)
                fitFitted(this);
            end
            
            this.updatePositions();
            
            % Now that we're ready, set visible on
            set(this.Handle, 'Visible', 'on');  
        end
        
        function config = get.Configuration(this)
            config = sftoolgui.FitFigureConfiguration(this.FitUUID);
            config.LegendOn = this.LegendOn;
            config.Grid = this.GridState;
            
            config.FittingConfig = sftoolgui.FittingConfiguration;
            config.FittingConfig.Visible = this.HFittingPanel.Visible;
            
            config.ResidualsConfig = sftoolgui.ResidualsConfiguration;
            config.ResidualsConfig.Visible = this.HResidualsPanel.Visible;
            
            config.ResultsConfig = sftoolgui.ResultsConfiguration;
            config.ResultsConfig.Visible = this.HResultsPanel.Visible;
            
            config.SurfaceConfig = sftoolgui.SurfaceConfiguration;
            config.SurfaceConfig.Visible = this.HSurfacePanel.Visible;
            config.SurfaceConfig.PredictionLevel = this.HSurfacePanel.PredictionLevel;
            
            config.ContourConfig = sftoolgui.SurfaceConfiguration;
            config.ContourConfig.Visible = this.HContourPanel.Visible;
        end
        
        function set.GridState (this, value)
            if ~(strcmpi(value, 'on') || strcmpi(value, 'off'))
                error(message('curvefit:FitFigure:InvalidGridStateValue'));
            end
            this.GridState = value;
            notify(this, 'GridStateChanged');
            notify(this, 'SessionChanged');
            iUpdateGrids(this);
        end
        
        function set.LegendOn (this, value)
            if ~islogical(value)
                error(message('curvefit:FitFigure:InvalidLegendOnValue'));
            end
            this.LegendOn = value;
            notify(this, 'LegendStateChanged');
            notify(this, 'SessionChanged');
            updateLegends(this);
        end
    end
    
    methods(Access = public)
        generateMCode(this, mcode);
        plotData(this, axesNeedPreparation);
        updateFitFigure(this, config);
    end
    
    methods(Hidden)
        % updateValidationData is public for testing purposes.
        updateValidationData(this);
        % isReady is public for testing purposes.
        tf = isReady(this);
    end
    
    methods(Access = protected)
        tf = isAssociatedFit(this, fit);
    end
    
    methods(Access = private)
        adjustmenu(this);
        calculateAndSetLimits(this);
        createToolbar(this);
        fitFigureDeleteFcn(this, src, event);
        fitFitted(this);
        fitNameUpdated(this);
        n = numberOfVisiblePlots(this);
        plotFitAndResids(this);
        setNoDataVisible(this);
        showAxisLimitsDialog(this, src, event);
        toggleExcludeMode(this, src, event);
        toggleGridState(this, src, event);
        toggleLegendState(this, src, event);
        toggleProperty(this, property);
        updateFittingData(this, axesNeedPreparation);
        updateFitTypeValues (this);
        updateGrids(this);
        updateInformationPanel(this);
        updateLegends(this);
        updatePlotControls(this, controls, selectedProperty);
        updatePositions(this);
        updateResultsArea(this);
    end
end

function this = iCreatePanels(this)
% iCreatePanels creates the plot panels, a panel container for them and a
% "no data" panel.

% Create the "container" for plot panels.
this.HPlotPanel = sftoolgui.PlotLayoutPanel(this.Handle);
this.HPlotPanel.Tag = 'sftoolPlotUIPanel';

% Create the plot panels
this.HSurfacePanel = sftoolgui.SurfacePanel(this, this.HPlotPanel.HUIPanel);
this.HSurfacePanel.Visible = 'on';
this.HResidualsPanel = sftoolgui.ResidualsPanel(this, this.HPlotPanel.HUIPanel);
this.HResidualsPanel.Visible = 'off';
this.HContourPanel = sftoolgui.ContourPanel(this, this.HPlotPanel.HUIPanel);
this.HContourPanel.Visible = 'off';

this.HPlotPanel.setPanels(this.HSurfacePanel, this.HResidualsPanel, this.HContourPanel);

% The order of the panels here determines the order the panels will show up
% in the menu and on the toolbar
this.PlotPanels = {this.HSurfacePanel, this.HResidualsPanel, ...
    this.HContourPanel};

% Create the "No Data" panel
this.HNoDataPanel = sftoolgui.MessagePanel(this.Handle, 'Select data to fit curves or surfaces.');
this.HNoDataPanel.Visible = 'on';

this.HFittingPanel = sftoolgui.EditFitPanel(this.Handle, ...
    this.HFitdev, this.FitUUID, this);
this.HFittingPanel.Tag = 'FittingUIPanel';

this.HResultsPanel = sftoolgui.InfoAndResultsPanel(this.Handle);
this.HResultsPanel.Tag = 'ResultsUIPanel';

this.HLimitSpinnersDialog = sftoolgui.LimitSpinnersDialog(this);
end

function iMakeFigureGood(this)
% iMakeFigureGood
figH = this.Handle;
hFitdev = this.HFitdev;
set(figH, 'Name', hFitdev.FitName);

createToolbar(this);
adjustmenu(this);

% Set rotate3d to be the default mode
r = rotate3d( figH );
hM = uigetmodemanager(figH);
set(hM,'DefaultUIMode','Exploration.Rotate3d');
sftoolgui.util.allowLegendDragInRotate3D( figH, r );

set(figH, 'ResizeFcn', @(src, evt) updatePositions( this ) );
end

function iUpdateGrids(fitFigure)
updateGrid(fitFigure.HResidualsPanel);
updateGrid(fitFigure.HSurfacePanel);
updateGrid(fitFigure.HContourPanel);
end
