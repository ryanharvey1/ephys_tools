function sfUpdateLegend(hAxes, state)
%SFUPDATELEGEND Update sftool legends

%   $Revision: 1.1.6.7 $  $Date: 2011/05/09 00:39:52 $
%   Copyright 2008-2010 The MathWorks, Inc.

% This method is called both to turn the legend on or off and
% to update the legend when plots changed. deleting/and
% recreating the legend seems to be the way to update because legends
% do not automatically update when new plots are added or removed.

hLegend = legend(hAxes);
legendInfo = {'Location', 'NorthEast'};
legendPos = [];
if ~isempty(hLegend)
    legendInfo = cfswitchyard( 'cfgetlegendinfo', hLegend);
    location = get(hLegend, 'Location');
    if strcmpi(location, 'none');
        legendPos = get(hLegend, 'Pos');
    end
    delete(legend(hAxes));
end
if state  % turn on the legend
    hLegend = iShowLegend( hAxes );
    
    if ~isempty(hLegend)
        set(hLegend, 'Interpreter', 'none');
        iFixContextMenu(hLegend);
        set(hLegend, legendInfo{:});
        % legendInfo (returned from cfgetlegendinfo) does not return a location
        % if it is 'none', which is the state if users manually move the
        % legend. Check location, if it is none, use the legendPos.
        if ~isempty(legendPos)
            set(hLegend, 'Location', 'none', 'Position', legendPos)
        end
    end
end
end

function hLegend = iShowLegend( hAxes )
% iShowLegend -- Show the legend for the given axes.
%
% This function is a wrapper around "legend( hAxes, 'show' )" with protection
% against some warnings and errors that are not relevant to SFTOOL.

% Suppress warning about empty plot
ws = warning( 'off', 'MATLAB:legend:PlotEmpty' );
warningCleanup = onCleanup( @() warning( ws ) );

try
    hLegend = legend( hAxes, 'show' );
catch anException
    % Might get an error about "Unsupported ARRAYFUN input type: handle"
    % We can ignore this and return an empty for the legend handle
    warning(message('curvefit:sftoolgui:sfUpdateLegend:ErrorShowingLegend', anException.identifier, anException.message));
    hLegend = [];
end

end

function iFixContextMenu(hLegend)
% The legend gets created with a context menu. However this context menu
% has some features that have a destructive affect on SFTOOL. In this
% function, we remove those features
cmh = get( hLegend, 'UIContextMenu' );
if isempty( cmh )
    % If there is no context menu, then there is no need to remove anything from
    % it.
    return
end
% The children (menu entries) of the context menu are hidden so we need
% to get around that
h = allchild( cmh );

% Our actions are based on tags of items that appear in the context menu so
% we need to get all of those tags.
tags = get( h, 'Tag' );

% Delete the entries that cause bad things to happen
TAGS_TO_DELETE = {'scribe:legend:mcode', 'scribe:legend:propedit', ...
    'scribe:legend:interpreter', 'scribe:legend:delete', ...
    'scribe:legend:refresh'};
tf = ismember( tags, TAGS_TO_DELETE );
delete( h(tf) );

% Need to delete the separator that comes before "Color..." as we have
% deleted the item ('Delete') that comes just before it
tf = ismember( tags, 'scribe:legend:color' );
set( h(tf), 'Separator', 'off' );
end
