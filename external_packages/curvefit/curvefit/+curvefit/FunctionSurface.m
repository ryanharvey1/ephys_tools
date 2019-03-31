classdef FunctionSurface < handle
    %FUNCTIONSURFACE An HG surface representation of a surface fit object
    %
    %   A FUNCTIONSURFACE is wrapper around an HG surface. This HG surface is a
    %   plot of a surface fit object and will respond to changes in that fit
    %   object or to changes in the axes limits.
    %
    %   Example:
    %       % Fit a surface to data
    %       fo = fit( [Ne, OT], Gn, 'polynomial' );
    %
    %       % Create a function surface attached to the current axes (GCA)
    %       h = FunctionSurface;
    %
    %       % Setting the FitObject property produces a plot
    %       h.FitObject = fo;
    %
    %       % Changing the axes limits causes the surface plot to update
    %       set( gca, 'XLim', [683, 3045], 'Ylim', [0, 100] );
    %
    %       % Changing the Fit Object in the Function surface also causes the
    %       % plot to update
    %       fo = fit( [Ne, OT], Gn, 'GRIDDATA' );
    %       h.FitObject = fo;

    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.7.2.1 $    $Date: 2011/08/01 00:51:22 $ 

    properties(SetAccess = 'private', GetAccess = 'private', Dependent = true)
        % Parent -- This is the axes in which the surface lives. It is a
        % dependent property and is inferred from the HgSurface property
        Parent
    end
       
    properties(SetAccess = 'private', GetAccess = 'private')
        % AxesLimitListeners -- Listeners on axes limits which call the
        % redraw method (obj.Parent).
        AxesLimitListeners
        % AxesListeners -- Listeners on the axes (obj.Parent).
        AxesListeners
        % SurfaceListeners -- Listeners on the surface (obj.HgSurface).
        SurfaceListeners
        % HgSurface -- This is the HG surface object that we are wrapping up
        HgSurface
        % HgLower -- This is the HG representation of the lower bound 
        HgLower
        % HgUpper -- This is the HG representation of the upper bound 
        HgUpper
    end
    properties(SetAccess = 'public', GetAccess = 'public', Dependent = true)
        % DisplayName -- This is the name that will get used if the surface is
        % displayed in a legend.
        DisplayName
    end
    properties(SetAccess = 'public', GetAccess = 'public')
        % FitObject -- This is the surface fit object that this FunctionSurface
        % is a representation of.
        FitObject
        % PredictionBounds -- Option to turn 'on' or 'off' the plotting of
        % prediction bounds with the surface
        PredictionBounds = 'off';
        % PredictionBoundsOptions -- Options for computing the prediction
        % bounds.
        PredictionBoundsOptions = curvefit.PredictionIntervalOptions;
    end
    
    methods
        function obj = FunctionSurface( hAxes )
            % FUNCTIONSURFACE Create a function surface object
            %
            %   H = FUNCTIONSURFACE() is a function surface object attached to
            %   the current axes (GCA).
            %
            %   H = FUNCTIONSURFACE( anAxes ) is a function surface object attached
            %   to the given axes, anAxes.
            if ~nargin
                hAxes = gca;
            end
            
            % Create HG surfaces
             obj.HgSurface = curvefit.Surface( hAxes, ...
                'XData', [], 'YData', [], 'ZData', [], 'CData', [], ...
                'XLimInclude', 'off', 'YLimInclude', 'off', ...
                'Tag', 'curvefit.FunctionSurface' );            
             obj.HgLower = curvefit.Surface( hAxes, ...
                'XData', [], 'YData', [], 'ZData', [], 'CData', [], ...
                'XLimInclude', 'off', 'YLimInclude', 'off', ...
                'Tag', 'curvefit.FunctionSurface.Lower' );
             obj.HgUpper = curvefit.Surface( hAxes, ...
                'XData', [], 'YData', [], 'ZData', [], 'CData', [], ...
                'XLimInclude', 'off', 'YLimInclude', 'off', ...
                'Tag', 'curvefit.FunctionSurface.Upper' );

            % Make the FunctionSurface visible on the legend
            setLegendable( obj )
            
            % Create Listeners
            createAxesListeners( obj );
            createAxesLimitListeners( obj );
            createSurfaceListeners( obj )
        end
        
        function set.FitObject( obj, fo )
            obj.FitObject = fo;
            redraw( obj );
        end
        
        function hAxes = get.Parent( obj )
            hAxes = get( obj.HgSurface, 'Parent' );
        end
        
        function n = get.DisplayName( obj )
            n = get( obj.HgSurface, 'DisplayName' );
        end
        function set.DisplayName( obj, n )
            set( obj.HgSurface, 'DisplayName', n );
            set( obj.HgLower,   'DisplayName', sprintf( 'Pred bnds (%s)', n ) );
        end
        
        function set.PredictionBounds( obj, value )
            if ischar( value ) && ismember( value, {'on', 'off'} ),
                obj.PredictionBounds = value;
                redraw( obj );
            else
                error(message('curvefit:FunctionSurface:InvalidPredictionBounds'));
            end
        end
        
        function set.PredictionBoundsOptions( obj, value )
            if isa( value, 'curvefit.PredictionIntervalOptions' )
                obj.PredictionBoundsOptions = value;
                redraw( obj );
            else
                error(message('curvefit:FunctionSurface:PredictionBoundsOptions'));
            end
        end
        
        function delete( obj )
            % DELETE
            %   DELETE( H ) deletes the Function Surface H including the HG
            %   surfaces that it holds. Use DELETEBUTKEEPSURFACE if you want to
            %   keep the HG surfaces but delete the function surface
            %
            %   See also DELETEBUTKEEPSURFACE.
            
            % Delete the surface listeners so that we don't get into a "deletion
            % loop"
            obj.SurfaceListeners = {};
            
            % Delete the HG surfaces
            delete( obj.HgSurface );
            delete( obj.HgLower );
            delete( obj.HgUpper );
        end
        
        function hSurface = deleteButKeepSurface( obj )
            %DELETEBUTKEEPSURFACE Delete the FunctionSurface without deleting the HG Surfaces
            %
            %   S = DELETEBUTKEEPSURFACE( H ) deletes the Function Surface H
            %   without deleting the HG surfaces that it holds. Instead the HG
            %   surfaces are returned as the output argument S. This allows the
            %   FunctionSurface to be used as a utility for plotting surfaces,
            %   e.g.,
            %
            %       fo = fit( [Ne, OT], Gn, 'polynomial' );
            %       % Now plot and get a regular HG surface
            %       hFS = FunctionSurface( hAxes );
            %       hFS.FitObject = fo;
            %       hSurf = deleteButKeepSurface( hFS );
            %
            %   S will be a vector of one or more handles to HG objects. For
            %   example, if bounds are on then S will include handles to the
            %   surfaces representing the upper and lower bounds.
            
            % Get the HG Surfaces 
            if strcmpi( obj.PredictionBounds, 'on' )
                hSurface = [obj.HgSurface; obj.HgLower; obj.HgUpper];
                % Ensure the HG surfaces aren't deleted
                obj.HgSurface = [];
                obj.HgLower   = [];
                obj.HgUpper   = [];
            else
                hSurface = obj.HgSurface;
                % Ensure the main HG surface isn't deleted
                obj.HgSurface = [];
            end
            % Set the axis limits to respond to the surface limits
            set( hSurface, 'XLimInclude', 'on', 'YLimInclude', 'on' );
            
            % Delete the FunctionSurface
            delete( obj );
        end
        
        function setAxesLimitListenersEnabled( obj, state )
            % setAxesLimitListenersEnabled   Activate/deactivate the axes listeners that
            % redraw the FunctionSurface.
             for k=1:length(obj.AxesLimitListeners)
                 obj.AxesLimitListeners{k}.Enabled = state;
             end
        end
    end
    
    methods(Access = 'private')
        function redraw( obj, ~, ~ )
            % redraw   Redraw a function surface
            %
            % Note: This method can be used as callback, i.e., the second and third arguments
            % are SRC and EVT.
            
            if isempty( obj.FitObject )
                % Update the surface object with the new values
                set( obj.HgSurface, 'XData', [], 'YData', [], 'ZData', [], 'CData', [] );
                set( obj.HgLower,   'XData', [], 'YData', [], 'ZData', [], 'CData', [] );
                set( obj.HgUpper,   'XData', [], 'YData', [], 'ZData', [], 'CData', [] );
            else
                % Form a grid for the surface.
                xlim = double(get( obj.Parent, 'XLim' ));
                ylim = double(get( obj.Parent, 'YLim' ));
                [xi, yi] = meshgrid( ...
                    linspace( xlim(1), xlim(2), 49 ), ...
                    linspace( ylim(1), ylim(2), 51 ) );
                
                if strcmpi( obj.PredictionBounds, 'on' )
                    % Evaluate the fit object and bounds over the grid
                    [zi, lb, ub] = iPredictionBounds( obj, xi, yi );
                    
                    % Update the surface object with the new values
                    set( obj.HgSurface, 'XData', xi, 'YData', yi, 'ZData', zi, 'CData', zi );
                    if isempty( lb ) % && isempty( ub )
                        set( obj.HgLower,   'XData', [], 'YData', [], 'ZData', [], 'CData', [] );
                        set( obj.HgUpper,   'XData', [], 'YData', [], 'ZData', [], 'CData', [] );
                    else
                        set( obj.HgLower,   'XData', xi, 'YData', yi, 'ZData', lb, 'CData', lb );
                        set( obj.HgUpper,   'XData', xi, 'YData', yi, 'ZData', ub, 'CData', ub );
                    end
                    
                else % strcmpi( obj.PredictionBounds, 'off' )
                    % Evaluate the fit object over the grid
                    zi = feval( obj.FitObject, xi, yi );
                    
                    % Update the surface object with the new values
                    set( obj.HgSurface, 'XData', xi, 'YData', yi, 'ZData', zi, 'CData', zi );
                    set( obj.HgLower,   'XData', [], 'YData', [], 'ZData', [], 'CData', [] );
                    set( obj.HgUpper,   'XData', [], 'YData', [], 'ZData', [], 'CData', [] );
                end
            end
        end
        
        function createAxesListeners( obj ) 
            obj.AxesListeners = {
                curvefit.listener( obj.Parent, 'ObjectBeingDestroyed', @obj.deleteCallback )
                };
        end
        
        function createAxesLimitListeners( obj )
            obj.AxesLimitListeners = {
                curvefit.proplistener( obj.Parent, 'XLim', @obj.redraw )
                curvefit.proplistener( obj.Parent, 'YLim', @obj.redraw )
                 };
        end
        
        function createSurfaceListeners( obj )
            obj.SurfaceListeners = {
                curvefit.listener( obj.HgSurface, 'ObjectBeingDestroyed', @obj.deleteCallback )
                curvefit.listener( obj.HgLower,   'ObjectBeingDestroyed', @obj.deleteCallback )
                curvefit.listener( obj.HgUpper,   'ObjectBeingDestroyed', @obj.deleteCallback )
                };
        end
        
        function deleteCallback( obj, ~, ~ )
            delete( obj );
        end
        
        function setLegendable( obj )
            % setLegendable   Make the FunctionSurface visible on the legend.
            %
            %   setLegendable( theFunctionSurface ) sets properties of the FunctionSurface so
            %   that it is correctly displayed on the legend. The upper bound is never
            %   displayed on the legend. The main surface and the lower bound are sometimes
            %   displayed on the surface.
           
            % Always hide upper bound from legend
            curvefit.setLegendable( obj.HgUpper, false );
            
            % Hide main surface and lower bound from legend if 'HGUsingMATLABClasses'
            if ~curvefit.isHandlegraphics( obj.HgSurface )
                curvefit.setLegendable( obj.HgSurface, false );
                curvefit.setLegendable( obj.HgLower, false );
            end
        end
    end

end

function [zi, lb, ub] = iPredictionBounds( obj, xi, yi)
% Evaluate the fit object and bounds over the grid
try
    [ci, zi] = predint( obj.FitObject, [xi(:), yi(:)], ...
        obj.PredictionBoundsOptions );
    lb = reshape( ci(:,1), size( xi ) );
    ub = reshape( ci(:,2), size( xi ) );
    zi = reshape( zi,      size( xi ) );
catch ME
    % We are looking to catch errors in PREDINT caused by an inability of the
    % fit object to compute bounds
    if  strcmp( ME.identifier, 'curvefit:predint:cannotComputePredInts' ) || ...
        strcmp( ME.identifier, 'curvefit:predint:missingInfo' );
        zi = feval( obj.FitObject, xi, yi );
        lb = [];
        ub = [];
    else
        rethrow( ME );
    end
end
end

