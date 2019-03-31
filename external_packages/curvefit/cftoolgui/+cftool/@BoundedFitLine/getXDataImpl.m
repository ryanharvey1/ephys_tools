function varargout = getXDataImpl(hObj, ~)
% getXDataImpl -- Implementation of get method for XData property
%
%   getXDataImpl(hObj, storedValue) is the value of the XData property for the
%   BoundedFitLine object hObj. The storedValue is ignored.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/04/11 20:30:39 $

varargout{1} = hObj.calcXYData();

end
