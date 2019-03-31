function h_ = plot( obj, varargin )
% PLOT  Plot a surface fit object.
%
%   PLOT(FO) plots the surface fit object FO over the range of the current axes,
%   if any, or otherwise over the range stored in the fit.
%
%   PLOT(FO, [X, Y], Z) plots Z versus X and Y and plots FO over the range of X
%   and Y.
%
%   PLOT(FO, [X, Y], Z, 'Exclude', EXCLUDEDATA) plots the excluded data in a
%   different color. EXCLUDEDATA is a logical array where true represents an
%   outlier.
%
%   PLOT(FO, ..., 'Style', STYLE) selects which way to plot the surface fit
%   object FO. STYLE may be any of the following strings
%
%       'Surface'   Plots the fit object as a surface (default)
%       'PredFunc'  Surface with prediction bounds for function
%       'PredObs'   Surface with prediction bounds for new observation
%       'Residuals' Plot the residuals (fit is the plane Z=0)
%       'Contour'   Make a contour plot of the surface
%
%   PLOT(FO, ..., 'Level', LEVEL) sets the confidence level to be used in the
%   plot. LEVEL is a positive value less than 1 and has a default of 0.95 (for
%   95% confidence). This option only applies to the 'PredFunc' and
%   'PredObs' plot styles.
%
%   PLOT(FO, ..., 'XLim', XLIM) and
%   PLOT(FO, ..., 'YLim', YLIM) sets the limits of the axes used for the plot.
%   By default the axes limits are taken from the data, XY. If no data is given,
%   then the limits are taken from the surface fit object, FO.
%
%   H = PLOT(FO, ...) returns a vector of handles of the plotted objects.
%
%   H = PLOT(FO, ..., 'Parent', HAXES) plots the fit object FO in the axes with
%   handle HAXES rather than the current axes. In addition the range is taken
%   from these axes rather than from the fit or the data unless specified by the
%   'XLim' or 'YLim' parameters.
%
%   See also: LINE, SURFACE, CFIT/PLOT.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.8 $    $Date: 2011/04/11 16:09:30 $

[XY, Z, in, out, style, level, hParent] = iParseInputs( obj, varargin{:} );

switch lower( style )
    case 'surface'
        h1 = iPlotSurface( hParent, obj );
        h2 = iPlotData( hParent, XY, Z, in, out );
        
    case 'predfunc'
        h1 = iPlotPredFunc( hParent, obj, level );
        h2 = iPlotData( hParent, XY, Z, in, out );
        
    case 'predobs'
        h1 = iPlotPredObs( hParent, obj, level );
        h2 = iPlotData( hParent, XY, Z, in, out );
        
    case {'residual', 'residuals'}
        if isempty( XY ) || isempty( Z )
            error(message('curvefit:surffit:sfit:plot:ResidualNeedsData'));
        end
        residuals = Z - feval( obj, XY );
        h1 = [];
        iPlotResidualReference( hParent );
        h2 = iStemPlotData( hParent, XY, residuals, in, out );
        if ~ishold( hParent )
            set( hParent, 'ZLim', iZLimitsForResidualsPlot( Z, residuals ) );
        end
        
    case 'contour'
        h1 = iPlotContour( hParent, obj );
        h2 = iPlotData( hParent, XY, zeros( size( Z ) ), in, out );
        
    otherwise
        error(message('curvefit:surffit:sfit:plot:InvalidStyle', style));
end

% Return arguments
if nargout
    h_ = [h1(:); h2(:)];
end

grid( hParent, 'on' );
box(  hParent, 'on' );

end

% -- Plot Functions for Different Styles
function h = iPlotSurface( hParent, obj )
% Plot the surface
hFS = curvefit.FunctionSurface( hParent );
hFS.FitObject = obj;
h = deleteButKeepSurface( hFS );
end

function h = iPlotPredFunc( hParent, obj, level )
% Plot the surface
hFS = curvefit.FunctionSurface( hParent );
hFS.PredictionBounds = 'on';
hFS.PredictionBoundsOptions.Interval = 'Functional';
hFS.PredictionBoundsOptions.Level = level;
hFS.FitObject = obj;
h = deleteButKeepSurface( hFS );
h = h(1);
end

function h = iPlotPredObs( hParent, obj, level )
% Plot the surface
hFS = curvefit.FunctionSurface( hParent );
hFS.PredictionBounds = 'on';
hFS.PredictionBoundsOptions.Interval = 'Observation';
hFS.PredictionBoundsOptions.Level = level;
hFS.FitObject = obj;
h = deleteButKeepSurface( hFS );
h = h(1);
end

function h = iPlotResidualReference( hParent )
xlim = get( hParent, 'XLim' );
ylim = get( hParent, 'YLim' );
h = patch( xlim([1,1,2,2,1]), ylim([1,2,2,1,1]), zeros( 1, 5 ), ...
    'FaceAlpha', 0.2, 'FaceColor', [0.2, 0.2, 0.2] );
% Hide this referencce plane from legend
curvefit.setLegendable( h, false );
end

function h = iPlotContour( hParent, obj )

xlim = get( hParent, 'XLim' );
ylim = get( hParent, 'YLim' );

[xi, yi] = meshgrid( ...
    linspace( xlim(1), xlim(2), 49 ), ...
    linspace( ylim(1), ylim(2), 51 ) );
zi = feval( obj, xi, yi );

[~, h] = contourf( xi, yi, zi );

end

function h = iPlotData( hParent, XY, Z, in, out )
% Plot the data, if any
h = [];
if ~isempty( XY ) && ~isempty( Z )
    if any( in )
        h(end+1) = line( XY(in,1), XY(in,2), Z(in), ...
            'LineStyle', 'none', ...
            'Marker', 'o', ...
            'MarkerEdgeColor', 'w', ...
            'MarkerFaceColor', 'b', ...
            'Parent', hParent );
    end
    if any( out )
        h(end+1) = line( XY(out,1), XY(out,2), Z(out), 'Parent', hParent );
        iSetExclusionLineProperties( h(end) );
    end
end

end

function h = iStemPlotData( hParent, XY, Z, in, out )
% Plot the data, if any, using a stem plot
h = [];
if ~isempty( XY ) && ~isempty( Z )
    holdState = iHold( hParent, 'on' );
    if any( in )
        h(end+1) = stem3( hParent, XY(in,1), XY(in,2), Z(in), 'filled' );
    end
    if any( out )
        h(end+1) = stem3( hParent, XY(out,1), XY(out,2), Z(out), 'rx-' );
        % Set the line width on the markers but not on the stems
        hX = findobj( get( h(end), 'Children' ), 'Marker', 'x' );
        iSetExclusionLineProperties( hX );
    end
    iHold( hParent, holdState );
end
end

% -- Misc Helper Functions
function iSetExclusionLineProperties( h )
set( h, ...
    'LineStyle', 'none', ...
    'Marker', 'x', ...
    'MarkerEdgeColor', 'r', ...
    'MarkerSize', 9, ...
    'LineWidth', 2 );
end

function zlim = iZLimitsForResidualsPlot( zData, residuals  )
% Compute the limits of the Z-axis a residuals plot.

% Magnitude of z-data
magnitudeOfZData = max( abs( zData ) );

% Start by estimating the limits as the range of the data
zlim = [min( residuals ), max( residuals )];

% Don't let the z-axes get too small
minZLim = max( 1e-8*magnitudeOfZData, 1e-5 );
% ... on either side of the z==0 plane.
zlim(1) = min( zlim(1), -minZLim );
zlim(2) = max( zlim(2),  minZLim );
end

function oldState = iHold( hParent, newState )
% iHold -- A wrapper around HOLD that returns the old state.

% Get the old state
if ishold( hParent ),
    oldState = 'on';
else
    oldState = 'off';
end

% Set the new state
hold( hParent, newState );
end

% -- Parse Inputs
function [XY, Z, in, out, style, level, hParent] = iParseInputs( obj, varargin )

% Parse inputs
if nargin >= 3 && isnumeric( varargin{1} ) && isnumeric( varargin{2} )
    XY = varargin{1};
    Z  = varargin{2};
    varargin = varargin(3:end);
else
    XY = zeros( 0, 2 );
    Z  = zeros( 0, 1 );
end

if size( XY, 2 ) ~= 2
    error(message('curvefit:surffit:sfit:plot:InvalidXY'));
end
if size( Z, 2 ) ~= 1
    error(message('curvefit:surffit:sfit:plot:InvalidZ'));
end
if size( XY, 1 ) ~= size( Z, 1 ),
    error(message('curvefit:surffit:sfit:plot:MismatchedXYAndZ'));
end

p = inputParser;
p.FunctionName = 'plot';
p.addParamValue( 'Exclude', [], @(v) islogical( v ) && numel( v ) == numel( Z ) );
p.addParamValue( 'Style', 'Surface' );
p.addParamValue( 'Level', 0.95, @(v) isnumeric( v ) && isscalar( v ) && 0 < v && v < 1 );
p.addParamValue( 'XLim',  [], @(v) isnumeric( v ) && numel( v ) == 2 && v(1) < v(2) );
p.addParamValue( 'YLim',  [], @(v) isnumeric( v ) && numel( v ) == 2 && v(1) < v(2) );
p.addParamValue( 'Parent', [], @(h) ishghandle( h, 'axes' ) );
p.parse( varargin{:} );

exclude = p.Results.Exclude;
style   = p.Results.Style;
level   = p.Results.Level;
xlim    = p.Results.XLim;
ylim    = p.Results.YLim;
hParent = p.Results.Parent;

if isempty( exclude )
    in  = true(  size( Z ) );
    out = false( size( Z ) );
else
    in = ~exclude;
    out = exclude;
end

% If "parent" is using default, then create a new axes to plot in.
if ismember( 'Parent', p.UsingDefaults )
    hParent = newplot;
    
    if ~ishold( hParent )
        view( hParent, 3 );
        
        % Get the axes limits from the user options, from the data or from the
        % the surface object
        xlim = iGetRange( xlim, XY(:,1), obj.xlim );
        ylim = iGetRange( ylim, XY(:,2), obj.ylim );
    end
end

% If we have been given limits then we need to set the axes limits
% appropriately.
if ~isempty( xlim )
    set( hParent, 'XLim', xlim );
end
if ~isempty( ylim )
    set( hParent, 'YLim', ylim );
end

end

function limits = iGetRange( userLimits, data, objLimits )
% The userLimits specify the limits to use.  By default the axes limits are
% taken from the data. If no data is given, then the limits are taken from the
% surface fit object.
if ~isempty( userLimits )
    limits = userLimits;
elseif ~isempty( data )
    limits = [min( data ), max( data )];
else
    limits = objLimits;
end

end
