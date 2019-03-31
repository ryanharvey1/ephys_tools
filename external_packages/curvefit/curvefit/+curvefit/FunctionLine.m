classdef FunctionLine < handle
    % FunctionLine   An HG representation of a curve fit object
    
    %   Copyright 2010-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.2.2.2 $    $Date: 2011/08/01 00:51:21 $
    
    %% Properties
    properties(Constant, GetAccess = 'private')
        % GRANULARITY   The number of points to add the user provided XData
        % to generate the XData for plotting.
        GRANULARITY = 293;
    end
    properties(SetAccess = 'private', GetAccess = 'private')
        % AxesListeners   Listeners on the axes (Parent of the Main Line).
        AxesListeners
        % LineListeners  Listeners on the main line and the upper and
        % lower bounds.
        LineListeners
        % MainLine   This is the line object that we are wrapping up. We
        % will use the function to redraw it in response changes in the
        % view, e.g., in response to changes in the axes limits.
        MainLine
        % LowerLine   A line that shows lower bounds on predictions using
        % the function.
        LowerLine
        % UpperLine   A line that shows upper bounds on predictions using
        % the function.
        UpperLine
    end
    properties(SetAccess = 'public', GetAccess = 'public', Dependent)
        % DisplayName   This is the name that will get used if the surface is
        % displayed in a legend.
        DisplayName = [];
        % LineWidth   The width of the line
        LineWidth = [];
        % LineStyle   The style of the line
        LineStyle = [];
        % Color   The color of the line
        Color = [];
    end
    properties(SetAccess = 'public', GetAccess = 'public')
        % FitObject   This is the curve fit object that this FunctionLine
        % is a representation of.
        FitObject = [];
        % XData   These are "special points" added to the lines plotted.
        % They are usually the x-data from the fit.
        XData = [];
        % PredictionBounds   Option to turn 'on' or 'off' the plotting of
        % prediction bounds of the curve
        PredictionBounds = 'off';
        % PredictionBoundsOptions   Options for computing the prediction
        % bounds.
        PredictionBoundsOptions = curvefit.PredictionIntervalOptions;
    end
    
    %% Get and set methods
    methods
        function name = get.DisplayName( obj )
            % get.DisplayName   Get the DisplayName
            
            % The DisplayName is stored in the MainLine
            name = get( obj.MainLine, 'DisplayName' );
        end
        function set.DisplayName( obj, name )
            % set.DisplayName   Set the DisplayName
            
            % The DisplayName is stored in the MainLine
            set( obj.MainLine,  'DisplayName', name );
            
            % We also need to update the DisplayName of lower bound with a
            % string derived from the DisplayName.
            set( obj.LowerLine, 'DisplayName', sprintf( 'Pred bnds (%s)', name ) );
        end
        
        function linewidth = get.LineWidth( obj )
            linewidth = get( obj.MainLine, 'LineWidth' );
        end
        function set.LineWidth( obj, linewidth )
            set( obj.MainLine, 'LineWidth',  linewidth );
        end
        
        function linestyle = get.LineStyle( obj )
            linestyle = get( obj.MainLine, 'LineStyle' );
        end
        function set.LineStyle( obj, linestyle )
            set( obj.MainLine, 'LineStyle',  linestyle );
        end
        
        function color = get.Color( obj )
            % get.Color   Get the color of a FunctionLine
            
            % The Color is a dependent property of the MainLine
            color = get( obj.MainLine, 'Color' );
        end
        function set.Color( obj, color )
            % set.Color   Set the color of a FunctionLine
            
            % The Color is a dependent property of the MainLine
            set( obj.MainLine, 'Color', color );
            % When we set the color of the main line we also need to set
            % the color of the bounds.
            
            set( obj.LowerLine, 'Color', color );
            set( obj.UpperLine, 'Color', color );
        end
        
        function set.FitObject( obj, fo )
            % set.FitObject   Set the FitObject
            obj.FitObject = fo;
            % Setting the FitObject requires that the line be redrawn.
            redraw( obj );
        end
        
        function set.XData( obj, xdata )
            % set.XData   Set the XData
            
            % Ensure that the XData is a row.
            obj.XData = xdata(:).';
            % Setting the XData requires that the line be redrawn.
            redraw( obj );
        end
        
        function set.PredictionBounds( obj, value )
            if ischar( value ) && ismember( value, {'on', 'off'} )
                obj.PredictionBounds = value;
                redraw( obj );
            else
                error(message('curvefit:FunctionLine:InvalidPredictionBounds'));
            end
        end
        
        function set.PredictionBoundsOptions( obj, value )
            if isa( value, 'curvefit.PredictionIntervalOptions' )
                obj.PredictionBoundsOptions = value;
                redraw( obj );
            else
                error(message('curvefit:FunctionLine:PredictionBoundsOptions'));
            end
        end
    end
    
    %% Public Methods
    methods
        function obj = FunctionLine( hAxes )
            % FunctionLine   Construct a FunctionLine.
            %
            %   curvefit.FunctionLine( hAxes ) constructs a FunctionLine that uses the given
            %   axes, hAxes, to plot in.
            %
            %   curvefit.FunctionLine() plots lines in GCA.
            %
            %   See also curvefit.FunctionSurface
            if ~nargin
                hAxes = gca;
            end
            
            % Create lines
            obj.MainLine  = iCreateLine( hAxes, 'curvefit.FunctionLine.MainLine' );
            obj.LowerLine = iCreateLine( hAxes, 'curvefit.FunctionLine.LowerLine' );
            obj.UpperLine = iCreateLine( hAxes, 'curvefit.FunctionLine.UpperLine' );
            
            % The lower and upper bounds should be dashed.
            set( obj.LowerLine, 'LineStyle', '--' );
            set( obj.UpperLine, 'LineStyle', '--' );
            
            % Create Listeners
            createAxesListeners( obj );
            createLineListeners( obj );
            
            % Always hide the upper line from the legend.
            curvefit.setLegendable( obj.UpperLine, false );
        end
        
        function delete( obj )
            % DELETE   Delete a FunctionLine
            %   DELETE( H ) deletes the Function Line H including the HG
            %   lines that it holds.
            
            % Delete the line listeners so that we don't get into a
            % "deletion loop"
            obj.LineListeners = {};
            
            % Delete the HG surfaces
            delete( obj.MainLine );
            delete( obj.LowerLine );
            delete( obj.UpperLine );
        end
    end
    
    %% Private methods
    methods(Access = 'private')
        function redraw( obj, ~, ~ )
            % redraw   Redraw the FunctionLine based on the current
            % FitObject and axes limits.
            %
            %   REDRAW( OBJ, SRC, EVT )
            hAxes = get( obj.MainLine, 'Parent' );
            xlim = double(get( hAxes, 'XLim' ));
            
            % The XData to use for the lines is the union of the XData
            % provided by the user and some linearly spaced points.
            xdata =  sort( [
                linspace( xlim(1), xlim(2), obj.GRANULARITY ), ...
                obj.XData
                ] );
            
            warnState = warning('off', 'all');
            cleanupObj = onCleanup(@() warning(warnState));
            if isempty( obj.FitObject )
                ydata = [];
                lb = [];
                ub = [];
                
            elseif strcmpi( obj.PredictionBounds, 'on' )
                [ydata, lb, ub] = iPredictionBounds( obj.FitObject, xdata, ...
                    obj.PredictionBoundsOptions );
                
            else
                ydata = feval( obj.FitObject, xdata(:) );
                ydata = ydata(:).';
                lb = [];
                ub = [];
            end
            
            iSetXYData( obj.MainLine,  xdata, ydata );
            iSetXYData( obj.LowerLine, xdata, lb );
            iSetXYData( obj.UpperLine, xdata, ub );
        end
        
        function createLineListeners( obj )
            obj.LineListeners = {
                curvefit.listener( obj.MainLine,  'ObjectBeingDestroyed', @obj.deleteCallback )
                curvefit.listener( obj.LowerLine, 'ObjectBeingDestroyed', @obj.deleteCallback )
                curvefit.listener( obj.UpperLine, 'ObjectBeingDestroyed', @obj.deleteCallback )
                };
        end
        
        function createAxesListeners( obj )
            hAxes = get( obj.MainLine, 'Parent' );
            obj.AxesListeners = {
                curvefit.proplistener( hAxes, 'XLim', @obj.redraw )
                };
        end
        
        function deleteCallback( obj, ~, ~ )
            % deleteCallback   Callback function for deletion events
            %
            %   deleteCallback( obj, src, evt )
            delete( obj );
        end
    end
end

%% Internal Functions
function [yi, lb, ub] = iPredictionBounds( fitObject, xi, options )
% iPredictionBounds   Evaluate the fit object and bounds
try
    [ci, zi] = predint( fitObject, xi, ...
        options.Level, options.Interval, options.Simultaneous );
    
    lb = reshape( ci(:,1), size( xi ) );
    ub = reshape( ci(:,2), size( xi ) );
    yi = reshape( zi,      size( xi ) );
catch ME
    % We are looking to catch errors in PREDINT caused by an inability of
    % the fit object to compute bounds
    if  strcmp( ME.identifier, 'curvefit:predint:cannotComputePredInts' ) ...
            || strcmp( ME.identifier, 'curvefit:predint:missingInfo' );
        yi = feval( fitObject, xi );
        lb = [];
        ub = [];
    else
        rethrow( ME );
    end
end
end

function iSetXYData( h, xdata, ydata )
% iSetXYData   Set the X- and Y-data on a line
%
% The Y-data may be empty, in which case the X-data will also be set to
% empty.
if isempty( ydata )
    set( h, 'XData', [], 'YData', [] );
else
    ydata = curvefit.nanFromComplexElements( ydata );
    set( h, 'XData', xdata, 'YData', ydata );
end
end

function hLine = iCreateLine( parent, tag )
% iCreateLine   Create a line with the given parent and tag.
hLine = line( 'Parent', parent, 'Tag', tag, ...
    'XData', [], 'YData', [], 'ZData', [], ...
    'XLimInclude', 'off' );
end
