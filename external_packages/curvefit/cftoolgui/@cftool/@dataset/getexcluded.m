function evec = getexcluded(ds,outlier)
%GETEXCLUDED Get an exclusion vector for this dataset/outlier combination

%   $Revision: 1.2.2.4 $  $Date: 2007/10/15 22:40:53 $
%   Copyright 2001-2007 The MathWorks, Inc.
NONE = cfswitchyard( 'cfGetNoneString' );

if ~isequal(outlier, NONE ) && ~isempty(outlier)
   % For convenience, accept either an outlier name or a handle
   if ischar(outlier)
      outlier = find(getoutlierdb,'name',outlier);
   end
   evec = cfswitchyard('cfcreateexcludevector',ds,outlier);
else
   evec = false(length(ds.x),1);
end
