function fitFitted( this )
%fitFitted FitFigure callback to Fitdev's FitFitted event

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/08/29 08:22:16 $

updateResultsArea(this);
if isFitted(this.HFitdev)
    plotFitAndResids(this);
else
    return;
end
end

