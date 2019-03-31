classdef SurfacePanel < sftoolgui.CloseablePanel
    %SURFACEPANEL A panel used by SFTOOL for plotting data and fits
    %
    %   SURFACEPANEL (fitFigure, parent)
    
    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.16.2.1 $    $Date: 2011/07/18 00:30:57 $

    properties(SetAccess = 'private', GetAccess = 'public')
        HAxes ;
        HSurfacePlot ;
        HCurvePlot ;
        HFittingDataLine ;
        HValidationDataLine ;
        HFittingExclusionLine ;
        HFittingInclusionLine ;
        HFitdev;
        HFitFigure;
    end
    
    properties(SetAccess = 'public', GetAccess = 'public', Dependent)
        % PredictionLevel is the level of confidence of the prediction
        % bounds displayed for the curve or surface. To turn off prediction
        % bounds use PredictionLevel = 0.
        PredictionLevel = 0;
    end
    
    properties(SetAccess = 'private', GetAccess = 'private')
        % HAxesViewController stores a handle to the AxesSwitchController
        HAxesViewController;
        % PrivatePredictionLevel is the private property associated with
        % the public dependent PredictionLevel property.
        PrivatePredictionLevel = 0;
        % Listeners stores a handle to various listeners
        Listeners
    end
    
    properties (Constant);
        Icon = 'surfaceSMALL.png';
        % Description is used as the toolbar button tooltip
        Description = xlate('Main plot');
        % Name is used as the menu label
        Name = xlate('&Main Plot');
    end
    
    methods
        function this = SurfacePanel(fitFigure, parent)
            this = this@sftoolgui.CloseablePanel(parent);
            this.HFitFigure = fitFigure;
            this.HFitdev = fitFigure.HFitdev;
            set(this.HUIPanel, 'Tag', 'SurfaceUIPanel');
            
            % Set up axes
            this.HAxes = axes('Parent', this.HUIPanel, ...
                'Tag', 'sftool surface axes', ...
                'Box', 'on', ...
                'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
            % Turn off the interpreter for all labels
            set( get( this.HAxes, 'XLabel' ), 'Interpreter', 'none');
            set( get( this.HAxes, 'YLabel' ), 'Interpreter', 'none');
            set( get( this.HAxes, 'ZLabel' ), 'Interpreter', 'none');
            
            % Set up graphics (surface & line) for the fit
            this.HSurfacePlot = curvefit.FunctionSurface(this.HAxes);
            this.HCurvePlot = curvefit.FunctionLine(this.HAxes);
            this.HCurvePlot.Color = 'r';
            this.HCurvePlot.LineWidth = 1.5;
            
            % Set up line for fitting data
            this.HFittingDataLine = sftoolgui.util.lineForExclusion( ...
                this.HAxes, 'SurfaceFittingDataLine' );
            
            % Set up lines to show other data: included, excluded, validation
            this.HFittingInclusionLine = line('XData', [], 'YData', [], ...
                'ZData', [], 'Parent', this.HAxes, ...
                'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', ...
                'w', 'MarkerFaceColor', 'b', ...
                'Tag', 'SurfaceFittingInclusionDataLine');
            this.HFittingExclusionLine = line('XData', [], 'YData', [], ...
                'ZData', [], 'Parent', this.HAxes, ...
                'LineStyle', 'none', 'Marker', 'pentagram', 'MarkerEdgeColor', ...
                'w', 'MarkerSize', 14, 'MarkerFaceColor', 'r', ...
                'Tag', 'SurfaceFittingExclusionDataLine');
            this.HValidationDataLine = line('XData', [], 'YData', [], ...
                'ZData', [], 'Parent', this.HAxes, ...
                'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', ...
                'b', 'MarkerFaceColor', 'w', ...
                'Tag', 'SurfaceValidationDataLine');
            
            % Create an AxesViewController
            this.HAxesViewController = ...
                sftoolgui.AxesViewController(this.HAxes, ...
                this.HFitFigure.AxesViewModel, ...
                isCurveDataSpecified(this.HFitdev.FittingData));
            
            % Set the limits
            limitsChangedAction(this, this.HFitFigure.AxesViewModel);
            
            % Create listeners
            createListeners(this);
        end
        
        function level = get.PredictionLevel(this)
            % Return the value of the associated private property
            level = this.PrivatePredictionLevel;
        end
        
        function set.PredictionLevel(this, level)
            % In addition to setting the PrivatePredictionLevel property,
            % this method will set the PredictionBounds and
            % PredictionBoundsOptions.Level properties of HSurfacePlot and
            % HCurvePlot. It will then (re)plot the curve or surface.
            this.PrivatePredictionLevel = level;
            if (level == 0)
                this.HSurfacePlot.PredictionBounds = 'off';
                this.HCurvePlot.PredictionBounds = 'off';
            else
                level = level/100;
                this.HSurfacePlot.PredictionBounds = 'on';
                this.HCurvePlot.PredictionBounds = 'on';
                this.HSurfacePlot.PredictionBoundsOptions.Level = level;
                this.HCurvePlot.PredictionBoundsOptions.Level = level;
            end
            plotSurface(this);
        end
        
        function plotSurface(this)
            clearSurface(this);
            if isCurveDataSpecified(this.HFitdev.FittingData)
                this.HCurvePlot.FitObject = this.HFitdev.Fit;
            else
                this.HSurfacePlot.FitObject = this.HFitdev.Fit;
            end
            updateDisplayNames(this);
            updateGrid(this);
            updateLegend(this, this.HFitFigure.LegendOn);
        end
        
        function clearSurface(this)
            this.HSurfacePlot.FitObject = [];
            this.HCurvePlot.FitObject = [];
        end
        
        function updateGrid(this)
            set(this.HAxes, 'XGrid', this.HFitFigure.GridState);
            set(this.HAxes, 'YGrid', this.HFitFigure.GridState);
            set(this.HAxes, 'ZGrid', this.HFitFigure.GridState);
        end
        
        function updateDisplayNames(this)
            this.HSurfacePlot.DisplayName = this.HFitdev.FitName;
            this.HCurvePlot.DisplayName = this.HFitdev.FitName;
            sftoolgui.util.refreshLegend(this.HAxes, this.HFitFigure.LegendOn)
        end
        
        function updateLegend(this, state)
            sftoolgui.sfUpdateLegend(this.HAxes, state);
        end
        
        function this = prepareAxesForPreview(this)
            % prepareAxesForPreview sets labels, updates the grid and sets
            % the AxesViewController View2D property
            [fx, fy, fz] = getValues(this.HFitdev.FittingData);
            [fxName, fyName, fzName] = getNames(this.HFitdev.FittingData);
            [vx, vy, vz] = getValues(this.HFitdev.ValidationData);
            [vxName, vyName, vzName] = getNames(this.HFitdev.ValidationData);
            set( get( this.HAxes, 'XLabel' ), ...
                'String', getLabel(fx, fxName, vx, vxName, 'X'));
            set( get( this.HAxes, 'YLabel' ), ...
                'String', getLabel(fy, fyName, vy, vyName, 'Y'));
            set( get( this.HAxes, 'ZLabel' ), ...
                'String', getLabel(fz, fzName, vz, vzName, 'Z'));
            updateGrid(this);
            
            this.HAxesViewController.View2D = isCurveDataSpecified(this.HFitdev.FittingData);
        end
        
        function plotDataLine(this, dLine, values, displayName)
            % if curve data is specified set zdata to empty. 
            if isCurveDataSpecified( this.HFitdev.FittingData )
                zdata = [];
            else
                zdata = values{3};
            end                   
            set(dLine, 'XData', values{1}, 'YData', values{2}, ...
                'ZData', zdata, 'DisplayName', displayName);
        end
        
        function plotDataLineWithExclusions(this, values, exclude, ...
                displayName)
            plotDataLine(this, this.HFittingDataLine, values, '');
            includedValues = values;
            includedValues{1} = values{1}(~exclude);
            includedValues{2} = values{2}(~exclude);
            includedValues{3} = values{3}(~exclude);
            plotDataLine(this, this.HFittingInclusionLine, includedValues, displayName);
            excludedValues = values;
            excludedValues{1} = values{1}(exclude);
            excludedValues{2} = values{2}(exclude);
            excludedValues{3} = values{3}(exclude);
            plotDataLine(this, this.HFittingExclusionLine, excludedValues, ...
                sprintf('Excluded %s',  displayName));
        end
        
        function clearFittingDataLine(this)
            clearDataLine(this.HFittingDataLine);
            clearDataLine(this.HFittingExclusionLine);
            clearDataLine(this.HFittingInclusionLine);
        end
        
        function clearValidationDataLine(this)
            clearDataLine(this.HValidationDataLine);
        end
        
        function tf = canGenerateCodeForPlot(this)	 
             % canGenerateCodeForPlot   True if code can be generated	 
             %	 
             %    canGenerateCodeForPlot( this ) returns true if code can be	 
             %    generated for surface plots and false otherwise. Code can	 
             %    be generated for surface plots if the panel is visible.	 
 	 
             tf = strcmpi( this.Visible, 'on' );	 
         end
        
        function generateMCode( this, mcode )
            % generateMCode   Generate code for a Surface Panel
            %
            %    generateMCode( H, CODE ) generates code for the given
            %    surface panel, H, and adds it the code object CODE.
            theFitdev = this.HFitdev;
            if isCurveDataSpecified( this.HFitdev.FittingData )
                cg = sftoolgui.codegen.CurvePlotCodeGenerator();
                
                xValues = getValues( theFitdev.FittingData );
                cg.HasXData = ~isempty( xValues );
            else
                cg = sftoolgui.codegen.SurfacePlotCodeGenerator();
                [cg.View(1), cg.View(2)] = view( this.HAxes );
            end
            
            cg.DoPredictionBounds = strcmpi( this.HSurfacePlot.PredictionBounds, 'on' );
            cg.PredictionLevel    = this.HSurfacePlot.PredictionBoundsOptions.Level;
            cg.FitName            = theFitdev.FitName;
            cg.FittingDataName    = theFitdev.FittingData.Name;
            cg.ValidationDataName = theFitdev.ValidationData.Name;
            cg.HaveLegend         = ~isempty( legend( this.HAxes ) );
            cg.HaveExcludedData   = any( theFitdev.Exclusions );
            cg.HaveValidation     = isValidationDataValid( theFitdev );
            cg.GridState          = this.HFitFigure.GridState;
            
            generateMCode( cg, mcode );
        end % of generateMCode
    end
    
    methods(Access = 'private')
        
        function dimensionChangedAction(this)
            % dimensionChangedAction sets the AxesViewController View2D
            % property.
            this.HAxesViewController.View2D = isCurveDataSpecified(this.HFitdev.FittingData);
        end
        
        function limitsChangedAction(this, model)
            % limitsChangedAction sets the axes limits. MODEL is
            % either an AxesViewModel or an AxesViewEventData.
            if ~isa(model, 'sftoolgui.AxesViewModel')
                model = model.AxesViewModel;
            end
            
            xlim = model.XInputLimits;
            if isCurveDataSpecified(this.HFitdev.FittingData)
                ylim = model.ResponseLimits ;
                zlim = [-1 1];
            else
                ylim = model.YInputLimits;
                zlim = model.ResponseLimits ;
            end
            set(this.HAxes, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim);
        end
        
        function createListeners(this)
            % createListeners adds listeners for the AxesViewModel
            % LimitsChanged event and the Fitdev DimensionChanged event.
            this.Listeners = {
                event.listener(this.HFitFigure.AxesViewModel, 'LimitsChanged', @(s, e) this.limitsChangedAction(  e ) );
                event.listener(this.HFitdev, 'DimensionChanged',  @(s, e) this.dimensionChangedAction ());
                };
        end
    end
end % of classdef

function clearDataLine(dLine)
set(dLine, 'XData', [], 'YData', [], 'ZData', []);
end

% Use data to determine name. If neither fitting nor validation data is
% specified, use the default label. If fitting data is specified use its
% name. If only validation data is specified, use its name.
function label = getLabel(fData, fName, vData, vName, defaultLabel)
label = defaultLabel;
if isempty(fData) && isempty(vData)
    return;
end
if ~isempty(fData)
    label = fName;
    return;
end
if ~isempty(vData)
    label = vName;
    return;
end
end

