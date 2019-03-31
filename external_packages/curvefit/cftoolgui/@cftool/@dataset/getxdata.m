function x = getxdata(ds,outlier)
%GETXDATA Get the X data not excluded by the specified excluded set

%   $Revision: 1.2.2.3 $  $Date: 2007/10/15 22:40:54 $
%   Copyright 2001-2007 The MathWorks, Inc.

NONE = cfswitchyard( 'cfGetNoneString' );

if ~isequal(outlier, NONE ) && ~isempty(outlier)
   % For convenience, accept either an outlier name or a handle
   if ischar(outlier)
      outlier = find(getoutlierdb,'name',outlier);
   end
   evec = cfswitchyard('cfcreateexcludevector',ds,outlier);
   x = ds.x(~evec);
else
   x = ds.x;
end
