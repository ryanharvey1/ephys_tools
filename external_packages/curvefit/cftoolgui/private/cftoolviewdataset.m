function [imsource, X, Y, W] = cftoolviewdataset( ds, width, height )
% For use by CFTOOL

%   [I, X, Y, W] = CFTOOLVIEWDATASET( DS, WIDTH, HEIGHT )
%
%   I -- image data for view plot
%   X, Y, W -- vectors of x and y data and weights
%   
%   DS -- curve fitting data set object
%   WIFTH, HEIGHT -- size of the image to return

%   $Revision: 1.1.6.2 $  $Date: 2008/09/13 06:48:12 $
%   Copyright 2000-2008 The MathWorks, Inc.

hDS = handle(ds);
X = hDS.x;
Y = hDS.y;
W = hDS.weight;

imsource = cftoolimagefromplot( @nPlot, width, height );

    function hLine = nPlot( hAxes )
        set( hAxes, ...
            'XLim', cfAxisLimitsFromData( X ), ...
            'YLim', cfAxisLimitsFromData( Y ), ...
            'xtick', [],...
            'ytick', [], ...
            'visible', 'off', 'View', [0, 90] );
        hLine = line( X, Y, 'parent', hAxes, 'marker', '.', 'linestyle' , 'none' );
    end
end
