function cfSetAxesLimits(limitName, limits)
%CFSETAXESLIMITS Set the limits of the main axes in CFTOOL
%
%   CFSETAXESLIMITS('XLim', XLIM)
%   CFSETAXESLIMITS('YLim', YLIM)
%
%   If the axes are not in manual limit mode (cfIsInManualLimitMode)
%   then the limits in the axes are not set. However the "reset limits"
%   are always set to the given XLIM or YLIM (unless the given value is
%   empty). Thus when a user selects "Default Axes Limits" the axes
%   return to the correct limits.
%
%   This function will also add a small non-zero margin to the limits.

%   Copyright 2007-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $    $Date: 2008/04/06 19:13:55 $
    
LIMIT_PROPERTIES = {'XLim', 'YLim'};

% If the given limits are empty then do nothing
if isempty( limits )
    return
end

% Get the handle of the main axes
cffig = cfgetset( 'cffig' );
hAxes = findall( cffig, 'Type', 'axes', 'Tag', 'main' );

% Are we in "manual limit mode"?
isInManualLimitMode = cfIsInManualLimitMode;

% Ensure that we have finite limits
if all( isfinite( limits ) )
    % Limits are finite

    % Add a small and non-zero margin to the limits
    dx = diff( limits ) * 0.01 * [-1 1];
    if all( dx==0 )
        dx = [-1 1] * max( 1, eps( limits(1) ) );
    end
    limits = limits + dx;

elseif isfinite( limits(1) )
    % limts(2) is NaN or Inf
    limits(2) = limits(1) + max( 1, eps( limits(1) ) );

elseif isfinite( limits(2) )
    % limts(1) is NaN or Inf
    limits(1) = limits(2) - max( 1, eps( limits(2) ) );

else
    % Both limits are NaN or Inf
    limits = [0, 1];
end
    
% Get the current x and y limits
if isInManualLimitMode
    oldLimits = get( hAxes, LIMIT_PROPERTIES );
end

% Set the axes limits to given limits
% -- We do a "zoom out" here to ensure that the "zoom reset" gets the
% correct limits for both the x- and y-axis. This is safe to do because
% either CFTOOL is in control of the limits and the axes limits will be
% at the "zoom out" location, or the user is in control of the limits
% and they will be reset below.
zoom( cffig, 'out' );
set( hAxes, limitName, limits );

% Call the "reset" function,
zoom( hAxes, 'reset' );

% Return the axes limits to their previous values
if isInManualLimitMode
    set( hAxes, LIMIT_PROPERTIES, oldLimits );
end
