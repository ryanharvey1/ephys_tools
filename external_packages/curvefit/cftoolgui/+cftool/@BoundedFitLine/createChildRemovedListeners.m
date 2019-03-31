function createChildRemovedListeners(hObj, hAxes)
%createChildRemovedListeners

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/04/11 20:30:34 $ 

if ~isempty(hObj.ChildRemovedListener) && isequal(hObj.ChildRemovedListener.Source,hAxes.ChildContainer)
    return;
end
hObj.ChildRemovedListener = event.listener( hAxes.ChildContainer, 'ObjectChildRemoved', ...
    @(es,ed) hObj.MarkDirty('limits') );


end
