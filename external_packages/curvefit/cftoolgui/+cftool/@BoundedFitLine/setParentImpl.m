function varargout = setParentImpl(hObj, hPar)
% setParentImpl -- Implementation of set method for Parent property

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/04/11 20:30:45 $ 

hObj.ChildAddedListener = event.listener(ancestor(hPar,'axes'),'MarkedDirty',...
    @(es,ed) hObj.MarkDirty('limits'));

hObj.createChildRemovedListeners(ancestor(hPar,'axes'));

hObj.ChildContainerListener = event.listener(ancestor(hPar,'axes'),'MarkedClean',...
    @(es,ed) hObj.createChildRemovedListeners(es));

varargout{1} = hPar;

end
