function listener = proplistener(obj, propname, callback)
%PROPLISTENER  Listener object for property PropertyPostSet events 
%
%   L = PROPLISTENER(obj, propertyName, callback)
%
%   Example:
%   If the handle object h needs to redraw when some axes have their x-limits
%   changed the listener may be created as follows.
%
%       h.Listener = curvefit.proplistener( hAxes, 'XLim', @h.redraw )
%
%   Note:
%   This function only creates listeners for property post-set events.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $    $Date: 2011/05/09 00:39:20 $ 


if curvefit.isHandlegraphics( obj )
    % Convert double to handle
    obj = handle( obj );
    % Find property on the first object
    property = findprop( obj(1), propname );
    % Make listener
    listener = handle.listener( obj, property, 'PropertyPostSet', callback );
else
    property = findprop( obj(1), propname );
    listener = event.proplistener( obj, property, 'PostSet', callback );
end

end
