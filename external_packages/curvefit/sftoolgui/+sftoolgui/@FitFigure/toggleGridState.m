function toggleGridState(fitFigure, ~, ~)
%toggleGridState Toggle grid state
%
%   toggleGridState(fitFigure, SOURCE, EVENT) is the callback to Grid menu
%   item and toolbar button click.

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/08/29 08:22:25 $

toggleProperty(fitFigure, 'GridState');
end