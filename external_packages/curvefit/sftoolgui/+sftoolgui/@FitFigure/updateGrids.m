function updateGrids(fitFigure)
%updateGrids Update all FitFigure plot panels' grids  
%
%   updateGrids is called to update all FitFigure plot panels' grids

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/08/29 08:22:31 $

updateGrid(fitFigure.HResidualsPanel);
updateGrid(fitFigure.HSurfacePanel);
updateGrid(fitFigure.HContourPanel);
end