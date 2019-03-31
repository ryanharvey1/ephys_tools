function cfSurface = Surface(hAxes, varargin)
% Surface    Factory function for surfaces
%
%   curvefit.Surface( anAxes, ... ) is a surface whose parent is anAxes. 

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2011/03/22 18:29:35 $

% if axes is a handle graphics object,
if curvefit.isHandlegraphics(hAxes)
    % then, use surface
    cfSurface = surface( 'Parent', hAxes, varargin{:} );
else
    % else, use a primitive surface
    cfSurface = primitiveSurface( hAxes, varargin{:} );
end
