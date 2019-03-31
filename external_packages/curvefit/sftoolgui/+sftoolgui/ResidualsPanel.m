classdef ResidualsPanel < sftoolgui.CloseablePanel
    %RESIDUALSPANEL A panel used by SFTOOL for residuals
    %
    %   RESIDUALSPANEL (fitFigure, parent)
    
    %   Copyright 2008-2011 The MathWorks, Inc.
    %    $Revision: 1.1.6.19.2.2 $    $Date: 2011/07/18 00:30:56 $

    properties(SetAccess = 'private', GetAccess = 'public')
        HAxes ;
        HResidualsPlot ;
        HResidualsLineForExclude ;
        HValidationDataPlot ;
        HReferencePlot;
        HExclusionPlot;
        HFitFigure ;
    end
    
    properties(SetAccess = 'private', GetAccess = 'private')
        % HAxesViewController stores a handle to the AxesViewController
        HAxesViewController;
        % Listeners -- various listeners
        Listeners
    end
    
    properties (Constant);
        Icon = 'residualSMALL.png';
        % Description is used as the toolbar button tooltip
        Description = xlate('Residuals plot');
        % Name is used as the menu label
        Name = xlate('Resi&duals Plot');
    end
       
    methods
        function this = ResidualsPanel(fitFigure, parent)
            this = this@sftoolgui.CloseablePanel(parent);
            this.HFitFigure = fitFigure;
            % Default Resid Panel should be invisible
            
            set(this.HUIPanel, 'Tag', 'ResidualUIPanel');
            
            this.HAxes = axes('Parent', this.HUIPanel, ...
                'Tag', 'sftool residuals axes', ...
                'Box', 'on', ...
                'NextPlot', 'add', ...
                'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
            
            this.HResidualsPlot = stem3([], [], [], 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'b', 'Parent', this.HAxes, 'Tag', 'ResidualsPlot');
            curvefit.setLegendable( this.HResidualsPlot, false );
            
            this.HExclusionPlot = stem3([], [], [], ...
                'Parent', this.HAxes, 'Marker', 'pentagram',  ...
                'MarkerSize', 14, 'Tag', 'ResidualsExclusionPlot', ...
                'Color', 'r', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r');
            curvefit.setLegendable( this.HExclusionPlot, false );
            
            % Turn off the interpreter for all labels
            set( get( this.HAxes, 'XLabel' ), 'Interpreter', 'none');
            set( get( this.HAxes, 'YLabel' ), 'Interpreter', 'none');
            set( get( this.HAxes, 'ZLabel' ), 'Interpreter', 'none');
            
            % vertexpicker (used in excluding outliers) can't handle a
            % stemseries. This artificial line will be used to select
            % excluded points
            this.HResidualsLineForExclude = sftoolgui.util.lineForExclusion( ...
                this.HAxes, 'ResidualsLineForExclude' );
            
            this.HValidationDataPlot = stem3([], [], [], 'Parent', this.HAxes, 'Tag', 'ResidualsValidationPlot');
            curvefit.setLegendable( this.HValidationDataPlot, false );
            
            % Plot the reference plane
            [X,Y] = meshgrid(get(this.HAxes, 'XLim'), get(this.HAxes, 'YLim'));
            this.HReferencePlot = surf(this.HAxes, X, Y, zeros(size(X)), 'FaceAlpha', .2, 'FaceColor', [.2 .2 .2], 'HitTest', 'off');
            % Don't (ever) show the reference plane in the legend
            curvefit.setLegendable( this.HReferencePlot, false );
            
            % Create an AxesViewController
            this.HAxesViewController = ...
                sftoolgui.AxesViewController(this.HAxes, ...
                this.HFitFigure.AxesViewModel, ...
                isCurveDataSpecified(this.HFitFigure.HFitdev.FittingData));
            
            % Set the limits
            limitsChangedAction(this, this.HFitFigure.AxesViewModel);
            
            % Create Listeners
            createListeners(this);
        end
        
        function updateGrid(this)
            set(this.HAxes, 'XGrid', this.HFitFigure.GridState);
            set(this.HAxes, 'YGrid', this.HFitFigure.GridState);
            set(this.HAxes, 'ZGrid', this.HFitFigure.GridState);
        end
        
        function plotDataLineWithExclusions(this)
            updatePlot(this);
        end
        
        function prepareAxesForPreview(this)
            % prepareAxesForPreview sets labels
            
            % Set the labels on the axes.
            fittingData = this.HFitFigure.HFitdev.FittingData;
            [xName, yName, zName] = getNames(fittingData);
            
            iSetLabel( this.HAxes, 'XLabel', xName, 'X' );
            iSetLabel( this.HAxes, 'YLabel', yName, 'Y' );
            iSetLabel( this.HAxes, 'ZLabel', zName, 'Z' );
            
            % update the axes view
            this.HAxesViewController.View2D = isCurveDataSpecified(this.HFitFigure.HFitdev.FittingData);
        end
        
        function updateLegend(this, state)
            hFitdev = this.HFitFigure.HFitdev;
            fitObject = hFitdev.Fit;
            if state && ~isempty(fitObject)
                curvefit.setLegendable( this.HResidualsPlot, true );
                
                haveValidationData = isValidationDataValid( hFitdev );
                curvefit.setLegendable( this.HValidationDataPlot, haveValidationData );
                
                anyPointsExcluded = ~isempty(hFitdev.Exclusions) && any(hFitdev.Exclusions);
                curvefit.setLegendable( this.HExclusionPlot, anyPointsExcluded );
            end
            sftoolgui.sfUpdateLegend(this.HAxes, state && ~isempty(fitObject));
        end
        
        function plotResiduals(this)
            updatePlot(this);
            updateLegend(this, this.HFitFigure.LegendOn);
        end
        
        function updateDisplayNames(this)
            hFitdev = this.HFitFigure.HFitdev;
            set(this.HResidualsPlot, 'DisplayName', sprintf('%s - residuals', hFitdev.FitName));
            set(this.HExclusionPlot, 'DisplayName', sprintf('Excluded %s', hFitdev.FittingData.Name));
            set(this.HValidationDataPlot, 'DisplayName', sprintf('%s - validation residuals', hFitdev.FitName));
            sftoolgui.util.refreshLegend(this.HAxes, this.HFitFigure.LegendOn)
        end
        
        function plotValidationData(this)
            updatePlot(this);
            updateLegend(this, this.HFitFigure.LegendOn);
        end
        
        function tf = canGenerateCodeForPlot(this)
            % canGenerateCodeForPlot   True if code can be generated
            %            
            %    canGenerateCodeForPlot( this ) returns true if code can be
            %    generated for residual plots and false otherwise. Code can 
            %    be generated for residuals plots if the panel is visible.
            
            tf = strcmpi( this.Visible, 'on' );
        end
        
        function generateMCode( this, mcode )
            % GENERATEMCODE   Generate code for a Residuals Panel
            %
            %    GENERATEMCODE( H, CODE ) generates code for the given
            %    residuals panel, H, and adds it the code object CODE.
            hFitdev = this.HFitFigure.HFitdev;
            
            if isCurveDataSpecified( hFitdev.FittingData )
                cg = sftoolgui.codegen.CurveResidualPlotCodeGenerator();
                
                xValues = getValues( hFitdev.FittingData );
                cg.HasXData = ~isempty( xValues );
            else
                cg = sftoolgui.codegen.SurfaceResidualPlotCodeGenerator();
                [cg.View(1), cg.View(2)] = view( this.HAxes );
            end
            
            cg.FitName            = this.HFitFigure.HFitdev.FitName;
            cg.FittingDataName    = this.HFitFigure.HFitdev.FittingData.Name;
            cg.ValidationDataName = this.HFitFigure.HFitdev.ValidationData.Name;
            cg.HaveLegend         = ~isempty( legend( this.HAxes ) );
            cg.HaveExcludedData   = any( this.HFitFigure.HFitdev.Exclusions );
            cg.HaveValidation     = isValidationDataValid(this.HFitFigure.HFitdev);
            cg.GridState          = this.HFitFigure.GridState;
            generateMCode( cg, mcode );
        end
        
        function [resids, vResids] = getResiduals(this)
            % getResiduals gets residuals and validation data residuals for
            % either curve or surface data.
            hFitdev = this.HFitFigure.HFitdev;
            fitObject = hFitdev.Fit;
            
            resids = [];
            vResids = [];
            
            if ~isempty(fitObject)
                resids = iResiduals( fitObject, hFitdev.FittingData );
                if isValidationDataValid(hFitdev)
                    vResids = iResiduals( fitObject, hFitdev.ValidationData );
                end
            end
        end
    end
    
    methods(Access = 'private')
        
        function createListeners(this)
            this.Listeners = {
                event.listener(this.HFitFigure.AxesViewModel, 'LimitsChanged', @(s, e) this.limitsChangedAction(  e ) );
                curvefit.listener( this.HAxes, 'ObjectBeingDestroyed', @(s, e) delete( this ) )
                event.listener(this.HFitFigure.HFitdev, 'DimensionChanged',  @(s, e) this.dimensionChangedAction ());
                };
        end
        
        function dimensionChangedAction(this)
            % dimensionChangedAction updates the axesViewController's
            % View2D property.
            this.HAxesViewController.View2D = isCurveDataSpecified(this.HFitFigure.HFitdev.FittingData);
        end
        
        function limitsChangedAction(this, model)
            % limitsChangedAction sets the limits (even if there are no
            % residuals). MODEL is either an AxesViewModel or an
            % AxesViewEventData.
            
            if ~isa(model, 'sftoolgui.AxesViewModel')
                model = model.AxesViewModel;
            end
            
            xlim = model.XInputLimits;
            if isCurveDataSpecified(this.HFitFigure.HFitdev.FittingData)
                ylim = model.ResidualLimits ;
                zlim = [-1 1];
            else
                ylim = model.YInputLimits;
                zlim = model.ResidualLimits ;
            end
            
            set(this.HAxes, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim);
        end
        
        function updateReferencePlane(this, ~, ~)
            % updateReferencePlane(this, src, evt)
            
            hFitdev = this.HFitFigure.HFitdev;
            % For surfaces, we want to plot a reference plane over axes limits
            [X,Y] = meshgrid(get(this.HAxes, 'XLim'), get(this.HAxes, 'YLim'));
            % For curves, where we want the "plane" to look like a "line",
            % set Y to zeros.
            if isCurveDataSpecified(hFitdev.FittingData)
                Y = zeros(size(X));
            end
            set(this.HReferencePlot, 'XData', X, 'YData', Y, 'ZData', zeros(size(X)));
        end
        
        function updatePlot(this)
            if strcmpi(this.Visible, 'on')
                set(this.HResidualsPlot, 'XData', [], 'YData', [], 'ZData', []);
                set(this.HValidationDataPlot, 'XData', [], 'YData', [], 'ZData', []);
                set(this.HExclusionPlot, 'XData', [], 'YData', [], 'ZData', []);
                set(this.HResidualsLineForExclude, 'XData', [], 'YData', [], 'ZData', []);
                hFitdev = this.HFitFigure.HFitdev;
                fitObject = hFitdev.Fit;
                
                [resids, vResids] = getResiduals(this);
                % Nothing is plotting unless there is a valid fit.
                if ~isempty(fitObject)
                    if isCurveDataSpecified(hFitdev.FittingData)
                        plotCurveResiduals(this, resids, vResids)
                    else
                        plotSurfaceResiduals(this, resids, vResids)
                    end
                end
                updateDisplayNames(this);
                updateReferencePlane(this);
            end
        end
        
        function plotCurveResiduals(this, resids, vResids)
            % Plot curve residuals
            hFitdev = this.HFitFigure.HFitdev;
            x = getCurveValues(hFitdev.FittingData);
            exclude = hFitdev.Exclusions;
            
            xInclude = x;
            xInclude(exclude) = NaN;
            yInclude = resids;
            yInclude(exclude) = NaN;
            setResidualData(this, xInclude, yInclude, []);
            
            xExclude = x;
            xExclude(~exclude) = NaN;
            yExclude = resids;
            yExclude(~exclude) = NaN;
            setExclusionData(this, xExclude, yExclude, []);
            
            setResidualLineForExclude(this, x, resids, []);
            
            if ~isempty(vResids)
                vx = getCurveValues(hFitdev.ValidationData);
                setValidationResidualData(this, vx, vResids, []);
            end
        end
        
        function plotSurfaceResiduals(this, resids, vResids)
            % Plot surface residuals
            
            hFitdev = this.HFitFigure.HFitdev;
            [x, y] = getValues(hFitdev.FittingData);
            exclude = hFitdev.Exclusions;
            
            xInclude = x;
            xInclude(exclude) = NaN;
            yInclude = y;
            yInclude(exclude) = NaN;
            residInclude = resids;
            residInclude(exclude) = NaN;
            setResidualData(this, xInclude, yInclude, residInclude);
            
            xExclude = x;
            xExclude(~exclude) = NaN;
            yExclude = y;
            yExclude(~exclude) = NaN;
            residExclude = resids;
            residExclude(~exclude) = NaN;
            setExclusionData(this, xExclude, yExclude, residExclude);
            
            setResidualLineForExclude(this, x, y, resids);
            
            if ~isempty(vResids)
                [vx, vy] = getValues(hFitdev.ValidationData);
                setValidationResidualData(this, vx, vy, vResids);
            end
        end
        
        function setResidualData(this, x, y, z)
            set(this.HResidualsPlot, ...
                'XData', x, ...
                'YData', y, ....
                'ZData', z);
        end
        
        function setExclusionData(this, x, y, z)
            set(this.HExclusionPlot, ...
                'XData', x, ...
                'YData', y, ....
                'ZData', z);
        end
        
        function setResidualLineForExclude(this, x, y, z)
            set(this.HResidualsLineForExclude, ...
                'XData', x, ...
                'YData', y, ....
                'ZData', z);
        end
        
        function setValidationResidualData(this, x, y, z)
            set(this.HValidationDataPlot, ...
                'XData', x, ...
                'YData', y, ...
                'ZData', z);
        end
    end
end

function iSetLabel( hAxes, labelName, newValue, defaultValue )
% iSetLabel sets the axis label with newValue or the defaultValue if
% newValue is empty.
if isempty(newValue)
    newValue = defaultValue;
end
set( get( hAxes, labelName ), 'String', newValue );
end

function residuals = iResiduals (fitObject, data)
% iResiduals calculates residuals. This method assumes that fitObject is
% not empty.
[input, output] = iGetInputOutput( data );
residuals = output - fitObject( input );
end

function [input, output] = iGetInputOutput ( data )
% iGetInputOutput returns input and output for curve or non-curve data.
% This method assumes that either curve or surface data has been specified.
if isCurveDataSpecified(data)
    [input, output] = getCurveValues(data);
else
    [x, y, output] = getValues(data);
    input = [x, y];
end
end
