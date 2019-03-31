function imsource = cftoolmakepreview( xname, yname )
% For use by CFTOOL

%   I = CFTOOLMAKEPREVIEW( X, Y )
%
%   I -- image data for preview plot
%
%   X, Y -- names of variables in the MATLAB workspace
%     These may '(none)' to use default values

%   $Revision: 1.10.2.7 $  $Date: 2008/09/13 06:48:11 $
%   Copyright 2000-2008 The MathWorks, Inc.


NONE = cfGetNoneString();

% if only x is given, swap x and y
if isequal( yname, NONE )
    yname = xname;
    xname = NONE;
end

% Get y from the workspace
y = evalin( 'base', yname );

% Get x, either implicitly or from the workspace
if isequal( xname, NONE )
    x = 1:length( y );
else
    x = evalin( 'base', xname );
end

% Get the plot as an image
width = 250;
height = 250;

imsource = cftoolimagefromplot( @nPlot, width, height );

    function hLine = nPlot( hAxes )
        set( hAxes, ...
            'XLim', cfAxisLimitsFromData( x ), ...
            'YLim', cfAxisLimitsFromData( y ), ...
            'xtick', [],...
            'ytick', [], ...
            'visible', 'off', 'View', [0, 90] );
        hLine = line( x, y, 'parent', hAxes, 'marker', '.', 'linestyle' , 'none' );
    end
end
