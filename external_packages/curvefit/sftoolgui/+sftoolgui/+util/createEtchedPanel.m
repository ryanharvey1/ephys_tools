function hPanel = createEtchedPanel(hParent)
%createEtchedPanel Create a uipanel with etched border
%
%   createEtchedPanel(hParent) creates a uipanel that is correctly set up
%   with the 'etchedin' border type and pixel units.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/05/09 00:40:01 $

if graphicsversion(hParent, 'handlegraphics')
    % Width=1 borders are doubled for etched types
    bWidth = 1;
else
    % Borders are not doubled for etched types
    bWidth = 2;
end

hPanel = uipanel(...
    'Parent', hParent, ...
    'Units', 'pixels', ...
    'BorderType', 'etchedin', ...
    'BorderWidth', bWidth);
