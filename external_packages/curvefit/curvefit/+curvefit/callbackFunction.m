function cb = callbackFunction(fcn, varargin)
%CALLBACKFUNCTION    Function for use in callbacks
%
%   CB = CALLBACKFUNCTION(FCN, ...)  is a object that can be used as a callback.
%
%   Essentially, this function takes the variable part of the argument list and
%   passes them into the function, FCN, after the source and event.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/01/23 20:37:28 $ 

cb = @callback;

    function callback( src, evt )
        fcn( src, evt, varargin{:} );
    end
end
