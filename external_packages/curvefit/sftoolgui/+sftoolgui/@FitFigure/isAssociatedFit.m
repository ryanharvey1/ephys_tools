function tf = isAssociatedFit( this, fit )
%isAssociatedFit determines whether or not FIT is the associated fit.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/08/07 07:25:04 $
    tf = (fit == this.HFitdev);
end

