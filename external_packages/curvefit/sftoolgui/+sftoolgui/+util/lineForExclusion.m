function theLine = lineForExclusion( parent, tag )
% lineForExclusion   A line that can be used to select points for inclusion or
% exclusion.
%
%   theLine = lineForExclusion( parent, tag ) 
%
%       parent -- the axes that is to draw the line in
%       tag -- a tag for the line

%   Copyright 2008-2011 The MathWorks, Inc.
%    $Revision: 1.1.6.1 $    $Date: 2011/05/09 00:40:07 $

theLine = line(...
    'Parent', parent, ...
    'XData', [], 'YData', [], 'ZData', [], ...
    'LineStyle', 'none', 'Marker', '+', 'MarkerSize', 1, 'Color', 'b', ...
    'Tag', tag );

% Don't (ever) show this artificial line in the legend
curvefit.setLegendable( theLine, false );

end
