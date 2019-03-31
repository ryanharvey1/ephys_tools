function W = getPanelBorderWidth(hPanel)
%getPanelBorderWidth Calculate the width of the panel's border decoration
%
%   W = getPanelBorderWidth(PANEL) returns the visual width of the panel's
%   current border settings for a single edge.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/05/09 00:40:05 $

W = get(hPanel, 'BorderWidth');

Type = get(hPanel, 'BorderType');
if strncmp(Type, 'etched', 6)
    % 'etchedin' or 'etchedout'
    if graphicsversion(hPanel, 'handlegraphics')
        % Etched borders are double the width
        W = 2*W;
    else
        % Etched borders use the border width
    end
elseif strcmp(Type, 'none')
    % Width ignored for none border type
    W = 0;
else
    % All other borders use the border width
end
