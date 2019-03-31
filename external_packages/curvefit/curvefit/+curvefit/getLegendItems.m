function [theGraphics, theText] = getLegendItems(anAxes)
%GETLEGENDITEMS   The graphics items displayed by a legend.
%
%   GETLEGENDITEMS( anAxes ) is a vector of graphics items in anAxes that
%   are displayed on the legend for the that plot
%
%   [theGraphics, theText] = GETLEGENDITEMS(anAxes) returns the text
%   displayed on the legend as well as the graphics reflected on the
%   legend. The text will be a cell array of strings with one element for each
%   graphic on the legend.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/04/11 16:09:27 $

if graphicsversion( anAxes, 'handlegraphics' )
    [theGraphics, theText] = iGetLegendItemsUsingHandleGraphics( anAxes );
else
    [theGraphics, theText] = iGetLegendItemsUsingMATLABClasses( anAxes );
end

end

function [theGraphics, theText] = iGetLegendItemsUsingHandleGraphics( anAxes )
[~, ~, theGraphics, theText] = legend( anAxes );

% Ensure that the text is a cell array
if isempty( theText )
    theText = {};
end
end

function [theGraphics, theText] = iGetLegendItemsUsingMATLABClasses(anAxes)

% Without the DRAWNOW the legend object doesn't have the handles of the objects
% that it is displaying
drawnow( 'expose' )

hLegend = legend( anAxes );

if isempty( hLegend ) || isempty( hLegend.UserData )
    theGraphics = [];
    theText = {};
else
    theGraphics = hLegend.UserData.handles;
    theText = hLegend.UserData.lstrings;
end

end
