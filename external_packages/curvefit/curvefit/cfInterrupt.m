function stop = cfInterrupt(action, stop)
%CFINTERRUPT   Interrupt the fitting process
%
%   STOP = CFINTERRUPT( 'get' ) gets the current interrupt status
%   CFINTERRUPT( 'set', STOP ) sets the current interrupt status.
%
%   STOP = true implies that fitting should/will be stopped.
%   STOP = false implies that fitting is free to continue.

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $    $Date: 2010/10/08 16:36:47 $

persistent PERSISTENT_STOP
if isempty( PERSISTENT_STOP )
    PERSISTENT_STOP = false;
end

switch action
    case 'get'
        stop = PERSISTENT_STOP;
    case 'set'
        PERSISTENT_STOP = stop;
    otherwise
        error(message('curvefit:cfInterrupt:InvalidAction'));
end

end
