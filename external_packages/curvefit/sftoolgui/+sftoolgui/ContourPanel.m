classdef ContourPanel < sftoolgui.CloseablePanel
    %CONTOURPANEL A panel used by SFTOOL for contour-plotting fits
    %
    %   CONTOURPANEL(fitFigure, parent)
    
    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.18.2.1 $    $Date: 2011/07/18 00:30:49 $
    
    properties(SetAccess = 'private', GetAccess = 'public')
        HAxes ;
        HContours ;
        HFittingDataLine ;
        HValidationDataLine ;
        HFittingExclusionLine ;
        HFitdev;
        HFitFigure;
    end
    
    properties (Constant);
        Icon = 'contourSMALL.png';
        % Description is used as the toolbar button tooltip
        Description = xlate('Contour plot');
        % Name is used as the menu label
        Name = xlate('&Contour Plot');
        NoCurveContoursMessage = message('curvefit:sftoolgui:ContourPanel:NoCurveContour');
    end
    
    properties (Constant, GetAccess = 'private');
        DefaultViewAngle = [0, 90];
    end
    
    properties(SetAccess = 'private', GetAccess = 'private')
        % UIPanel that contains the contour axes
        HImplementationPanel
        
        % MessagePanel object that displays when fit is 2D
        HMessagePanel
        
        % LimitsChangedListener listens for the LimitsChanged Event
        LimitsChangedListener
    end
    
    methods
        function this = ContourPanel( fitFigure, parent )
            this = this@sftoolgui.CloseablePanel(parent);
            set(this.HUIPanel, 'Tag', 'ContourUIPanel', 'BorderType', 'none');
            this.HFitFigure = fitFigure;
            this.HFitdev = fitFigure.HFitdev;
            this.HImplementationPanel = sftoolgui.util.createEtchedPanel(this.HUIPanel);
            set(this.HImplementationPanel, 'Tag', 'ContourImplementationUIPanel');
            
            % Create the contour and Axes
            createContourAndAxes(this);
            
            % Create the lines
            createTheLines(this);
            
            % Save this panel in the uipanel's appdata
            setappdata(this.HUIPanel, 'SurfaceFittingToolPlotPanel', this);
            
            % Create the MessagePanel
            this.HMessagePanel = sftoolgui.MessagePanel(this.HUIPanel, this.NoCurveContoursMessage.getString);
            this.HMessagePanel.Tag = 'ContourMessageUIPanel';
            
            this.LimitsChangedListener = event.listener(this.HFitFigure.AxesViewModel, ...
                'LimitsChanged', @(s, e) this.limitsChangedAction( e ) );
        end
        
        function plotSurface(this)
            iPlotContour( this );
            updateDisplayNames(this);
            updateGrid(this);
            updateLegend(this, this.HFitFigure.LegendOn);
        end
        
        function clearSurface(this)
            set( this.HContours, 'XData', [], 'YData', [], 'ZData', [] );
        end
        
        function updateGrid(this)
            set(this.HAxes, 'XGrid', this.HFitFigure.GridState);
            set(this.HAxes, 'YGrid', this.HFitFigure.GridState);
            set(this.HAxes, 'ZGrid', this.HFitFigure.GridState);
        end
        
        function updateDisplayNames(this)
            set( this.HContours, 'DisplayName', this.HFitdev.FitName );
        end
        
        function updateLegend(this, state)
            sftoolgui.sfUpdateLegend(this.HAxes, state);
        end
        
        % Use "plotting" (rather than real data) to get limits.
        % Plotting data contains zeros for undefined data values.
        function this = prepareAxesForPreview(this)
            updateVisibility(this);
            fData = this.HFitdev.FittingData;
            vData = this.HFitdev.ValidationData;
            
            [fx, fy] = getValues(fData);
            [fxName, fyName] = getNames(fData);
            [vx, vy] = getValues(vData);
            [vxName, vyName] = getNames(vData);
            set( get( this.HAxes, 'XLabel' ), ...
                'String', getLabel(fx, fxName, vx, vxName, 'X'));
            set( get( this.HAxes, 'YLabel' ), ...
                'String', getLabel(fy, fyName, vy, vyName, 'Y'));
            updateGrid(this);
        end
        
        function plotDataLine(~, dLine, values, displayName)
            %plotDataLine(this, dLine, values, displayName)
            % set z values to zero so that all points will show up on the
            % Contour plot.
            values{3} = zeros(size(values{3}));
            iPlotDataLine(dLine, values, displayName);
        end
        
        function plotDataLineWithExclusions(this, values, exclude, displayName)
            % set z values to zero so that all points will show up on the
            % Contour plot.
            values{3} = zeros(size(values{3}));
            
            includedValues{1} = values{1};
            includedValues{1}(exclude) = NaN;
            
            includedValues{2} = values{2};
            includedValues{2}(exclude) = NaN;
            
            includedValues{3} = values{3};
            includedValues{3}(exclude) = NaN;
            
            iPlotDataLine(this.HFittingDataLine, includedValues, displayName );
            
            excludedValues{1} = values{1};
            excludedValues{1}(~exclude) = NaN;
            
            excludedValues{2} = values{2};
            excludedValues{2}(~exclude) = NaN;
            
            excludedValues{3} = values{3};
            excludedValues{3}(~exclude) = NaN;
            
            iPlotDataLine(this.HFittingExclusionLine, excludedValues, sprintf('Excluded %s', displayName));
            iPlotContour(this);
        end
        
        function clearFittingDataLine(this)
            clearDataLine(this.HFittingDataLine);
            clearDataLine(this.HFittingExclusionLine);
        end
        
        function clearValidationDataLine(this)
            clearDataLine(this.HValidationDataLine);
            
        end
        
        function tf = canGenerateCodeForPlot(this)
            % CANGENERATECODEFORPLOT   True if code can be generated
            %
            %    CANGENERATECODEFORPLOT( this ) returns true if code can be
            %    generated for contour plots and false otherwise. Code can
            %    be generated for contour plots if the panel is visible and
            %    curve data is not specified.
            
            tf = strcmpi( this.Visible, 'on' ) && ...
                ~isCurveDataSpecified( this.HFitdev.FittingData );
        end
        
        function generateMCode( this, mcode )
            % GENERATEMCODE   Generate code for a Contour Panel
            %
            %    GENERATEMCODE( H, CODE ) generates code for the given
            %    contour panel, H, and adds it to the code object CODE.
            %
            %    This method should not be called when there is curve data
            %    selected. There is no check for curve data selected in
            %    in this method as this method should be called only if
            %    canGenerateCodeForPlot(this) is true.
            %
            %    See also canGenerateCodeForPlot
            
            cg = sftoolgui.codegen.ContourPlotCodeGenerator();
            cg.FitName            = this.HFitdev.FitName;
            cg.FittingDataName    = this.HFitdev.FittingData.Name;
            cg.ValidationDataName = this.HFitdev.ValidationData.Name;
            cg.HaveLegend         = ~isempty( legend( this.HAxes ) );
            cg.HaveExcludedData   = any( this.HFitdev.Exclusions );
            cg.HaveValidation     = isValidationDataValid(this.HFitdev);
            cg.GridState          = this.HFitFigure.GridState;
            generateMCode( cg, mcode );
        end % of generateMCode
    end
    
    methods(Access = protected)
        function layoutPanel(this)
            % Update the positions of objects within the panel
            
            this.layoutPanel@sftoolgui.CloseablePanel()
            
            % Set the Message panel and the Contour panel to match the
            % height and width of the Container panel.
            position = this.InnerPosition;
            
            this.HMessagePanel.Position = position;
            set(this.HImplementationPanel, 'Position', ...
                sftoolgui.util.adjustControlPosition(this.HImplementationPanel, position));
        end
    end
    
    methods(Access = private)
        function limitsChangedAction(this, evt)
            % limitsChangedAction updates the axes limits if curve data is
            % not specified.
            if ~isCurveDataSpecified(this.HFitdev.FittingData)
                xlim = evt.AxesViewModel.XInputLimits;
                ylim = evt.AxesViewModel.YInputLimits;
                zlim = evt.AxesViewModel.ResponseLimits ;
                % Contour plots must have 0 in the Z range
                zlim(1) = min(zlim(1), -1);
                zlim(2) = max(zlim(2), 1);
                set(this.HAxes, 'XLim', xlim, 'YLim', ylim, 'ZLim', zlim);
            end
        end
        
        function createContourAndAxes (this)
            this.HAxes = axes('Parent', this.HImplementationPanel );
            
            % Use the contour function to create a contour group.
            [~, this.HContours] = contour( this.HAxes, [] );
            set( this.HContours, ...
                'HitTest', 'off', ...
                'Fill', 'on', ...
                'LineColor', 'k', ...
                'Tag', 'sftoolgui.ContourPlot');
            curvefit.setLegendable( this.HContours, false );
            
            % Adding the contour group to axes changes some settings so we set
            % the things we want AFTER creating the contour group
            set( this.HAxes, ...
                'Tag', 'sftool contour axes', ...
                'Box', 'on', ...
                'XLim', [0 1], 'YLim', [0 1], 'ZLim', [-1 1], ...
                'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
                'View', this.DefaultViewAngle);
            % Don't allow Rotate3D on the contour plot
            behavior = hggetbehavior( this.HAxes, 'Rotate3d' );
            set( behavior, 'Enable', false );
            
            % Turn off the interpreter for x and y labels. (Z label is not
            % shown in contour plots).
            set( get( this.HAxes, 'XLabel' ), 'Interpreter', 'none');
            set( get( this.HAxes, 'YLabel' ), 'Interpreter', 'none');
            
            % Set the contour group to the bottom of the stack, to ensure
            % that the points are never hidden
            uistack( double( this.HContours ), 'bottom');
        end
        
        function createTheLines(this)
            this.HFittingDataLine = line( ...
                'Parent', this.HAxes, ...
                'XData', [], 'YData', [], 'ZData', [], ...
                'LineStyle', 'none', ...
                'Marker', 'o', ...
                'MarkerSize', 6, ...
                'MarkerFaceColor', 'b', ...
                'MarkerEdgeColor', 'w', ...
                'Tag', 'SurfaceFittingDataLine') ;
            
            this.HFittingExclusionLine = line('XData', [], 'YData', [], ...
                'ZData', [], 'Parent', this.HAxes, ...
                'LineStyle', 'none', 'Marker', 'pentagram', 'MarkerEdgeColor', ...
                'w', 'MarkerSize', 14,'MarkerFaceColor', 'r', ...
                'Tag', 'SurfaceFittingExclusionDataLine');
            
            this.HValidationDataLine = line( ...
                'Parent', this.HAxes, ...
                'XData', [], 'YData', [], 'ZData', [], ...
                'LineStyle', 'none', ...
                'Marker', 'o', ...
                'MarkerEdgeColor', 'b', ...
                'MarkerFaceColor', 'w', ...
                'Tag', 'SurfaceValidationDataLine' );
        end
        
        function updateVisibility(this)
            % updateVisibility sets the Visible property of the Message and
            % Contour panel depending on whether or not curve data is
            % specified.
            tf = isCurveDataSpecified(this.HFitdev.FittingData);
            this.HMessagePanel.Visible = iBooleanToOnOff( tf );
            set(this.HImplementationPanel, 'Visible', iBooleanToOnOff( ~tf ));
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

function iPlotContour( this )
if ~isCurveDataSpecified(this.HFitdev.FittingData)
    fitObject = this.HFitdev.Fit;
    if isempty( fitObject )
        set( this.HContours, 'XData', [], 'YData', [], 'ZData', [] );
    else
        xlim = get( this.HAxes, 'XLim' );
        ylim = get( this.HAxes, 'YLim' );
        
        [xi, yi] = meshgrid( ...
            linspace( xlim(1), xlim(2), 49 ), ...
            linspace( ylim(1), ylim(2), 51 ) );
        zi = fitObject( xi, yi );
        
        set( this.HContours, 'XData', xi, 'YData', yi, 'ZData', zi );
    end
end
end

function iPlotDataLine(dLine, values, displayName)
set(dLine, 'XData', values{1}, 'YData', values{2}, ...
    'ZData', values{3}, 'DisplayName', displayName);
end

function onOff = iBooleanToOnOff(tf)
% iBooleanToOnOff returns 'on' if tf is true and 'off' otherwise.
if tf
    onOff = 'on';
else
    onOff = 'off';
end
end

