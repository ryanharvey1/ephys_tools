function [imsource, hFigure] = cftoolimagefromplot( plotfun, width, height )
% CFTOOLIMAGEFROMPLOT For use by CFTOOL

%   IMSOURCE = CFTOOLIMAGEFROMPLOT( PLOTFUN, WIDTH, HEIGHT )
%
%   PLOTFUN should be a function that takes a handle to Axes and makes a plot in
%   those axes. It should return handles to anything plotted so that they can be
%   deleted from the axes after the image has been generated.

%   $Revision: 1.1.6.5 $  $Date: 2010/07/19 12:49:03 $
%   Copyright 2000-2010 The MathWorks, Inc.

hFigure = i_Figure( width, height );
hAxes = i_Axes( hFigure );

% If data has a complex part, it will spit a warning to the command line, so
% turn off warnings before plotting.
warnstate = warning( 'off', 'all' );
warningCleanup = onCleanup( @() warning( warnstate ) );

% Ask the caller to plot the data
hPlot = plotfun( hAxes );
% ... and then capture an image of this plot.
im = i_GetImage( hFigure );

% give the image a black edge
im(1,:,:)   = 0; 
im(end,:,:) = 0; 
im(:,1,:)   = 0; 
im(:,end,:) = 0;

% Convert to Java MemoryImageSource
imsource = im2mis( im );

% Delete the line so the axes are clear for the next plot
delete( hPlot );

end

function hFigure = i_Figure( width, height )
TAG = sprintf( 'helper figure for %s', mfilename );
% Look for a figure that has the correct tag
hFigure = findall( 0, 'type', 'figure', 'tag', TAG );
if isempty( hFigure )
    % No figures found
    % Create a new one
    hFigure = figure(...
        'Units', 'pixels',...
        'HandleVisibility', 'off', ...
        'IntegerHandle', 'off', ...
        'Visible', 'off',...
        'PaperPositionMode', 'auto', ...
        'Color', 'w', ...
        'Tag', TAG, ...
        'WindowStyle', 'Normal', ...
        'DockControls', 'off' );
elseif numel( hFigure ) > 1
    % Found more than one figure.
    % Use the first one
    hFigure = hFigure(1);
end
set( hFigure, 'Position', [1, 1, width, height] );

end

function hAxes = i_Axes( hFigure )
hAxes = i_FindAxes( hFigure );

if isempty( hAxes )
    % No axes found, create a new one
    hAxes = axes( 'Parent', hFigure, 'Tag', i_AxesTag() );
elseif numel( hAxes ) > 1
    % Found more than one axes, use the first one
    hAxes = hAxes(1);
end
% Make sure axes take up the whole figure
set( hAxes, 'Position', [.05, .05, .9, .9] );

end

function hAxes = i_FindAxes( hFigure )
% Look for the axes that have the correct tag
hAxes = findall( hFigure, 'type', 'axes', 'tag', i_AxesTag() );
end

function tag = i_AxesTag()
tag = sprintf( 'helper axes for %s', mfilename );
end

function im = i_GetImage( hFigure )
% i_GetImage -- Get an image from the given figure in way that works either with
% or without HG using MATLAB classes

if feature( 'HGUsingMATLABClasses' )
    im = print( hFigure, '-RGBImage', '-opengl', '-r0' );
else
    hAxes = i_FindAxes( hFigure );
    im = hardcopy( hAxes, '-dzbuffer', '-r0' );
end
end
