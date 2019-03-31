function varargout = getBoundLinesImpl(hObj, ~)
% getBoundLinesImpl -- Implementation of get method for BoundLines property
%
%   getBoundLinesImpl(hObj, storedValue)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/04/11 20:30:37 $ 

varargout{1} = [hObj.BoundedLineLowerHandle, hObj.BoundedLineUpperHandle];

end