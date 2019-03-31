function toggleLegendState(fitFigure, ~, ~)
%toggleLegendState Toggle legend state
%
%   toggleLegendState(fitFigure, SOURCE, EVENT) is the callback to the
%   Legend menu item and a toolbar button click.

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/08/29 08:22:26 $

toggleProperty(fitFigure, 'LegendOn');
end