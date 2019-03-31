function this = fitNameUpdated( this )
%fitNameUpdated FitFigure callback to Fitdev's FitNameUpdated event

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/08/29 08:22:17 $

set(this.Handle, 'Name', this.HFitdev.FitName);
updateDisplayNames(this.HResidualsPanel);
updateDisplayNames(this.HSurfacePanel);
updateDisplayNames(this.HContourPanel);
end