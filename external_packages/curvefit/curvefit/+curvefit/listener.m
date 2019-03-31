function listener = listener(obj, eventname, callback)
%LISTENER Create a listener on an event
%
%   L = LISTENER(obj, eventName, callback) creates a listener on a graphics
%   object.  The eventName should always be for the contemporary version of
%   the event, and will be mapped to legacy events automatically.
%
%   Examples: 
%
%   (1) If the handle object h holds a surface and should be deleted when
%   the surface is, then a listener can be set up like this:
%
%       h.Listener = curvefit.listener( h.Surface, 'ObjectBeingDestroyed', ...
%           @h.deleteCallback ) 
%
%   (2) To listen to resize events on a uipanel, a listener can be created
%   on 'SizeChange', no matter which version of the panel is used:
%       
%       h.Listener = curvefit.listener( uipanel, 'SizeChange', @h.resize ) 

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $    $Date: 2011/05/09 00:39:19 $ 


if curvefit.isHandlegraphics( obj )
    eventname = iTranslateEventName(eventname);
    listener = handle.listener( obj, eventname, callback );
else
    listener = event.listener( obj, eventname, callback );
end
end

function name = iTranslateEventName(name)
    % Map new event names to old ones
    switch name
        case 'SizeChange'
            name = 'ResizeEvent';
        otherwise
            % No name change is required
    end
end
