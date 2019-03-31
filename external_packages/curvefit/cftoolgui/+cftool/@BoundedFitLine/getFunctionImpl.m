function varargout = getFunctionImpl(hObj, ~)
% getFunctionImpl -- Implementation of get method for Function property
%
%   getFunctionImpl(hObj, storedValue) is the value of the Function property
%   for the BoundedFitLine object hObj. The storedValue is ignored.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/04/11 20:30:38 $

% The Function is a dependent property and corresponds to the "fit" property of
% the cftool.fit object held by this bounded line object.

if isempty( hObj.Fit )
    varargout{1} = [];
else
    varargout{1} = hObj.Fit.fit;
end

end
