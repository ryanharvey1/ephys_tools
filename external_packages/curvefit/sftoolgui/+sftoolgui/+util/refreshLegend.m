function refreshLegend( hAxes, legendOn )
% refreshLegend  refreshes the legend for the given axes if legendOn
%
%   refreshLegend( hAxes, legendOn )
%
%       hAxes -- the axes for the legend
%       legendOn -- legend state

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/05/09 00:40:08 $

% Refresh the legend if it is on and axes is not handlegraphics.
if legendOn && ~curvefit.isHandlegraphics(hAxes)
    sftoolgui.sfUpdateLegend(hAxes, true);
end

end
