function zoomModeCallback(fitFigure, src, ~, inout)
%zoomModeCallback Zoom mode callback
%
%   zoomModeCallback(fitFigure, SRC, EVENT) is the callback to the zoom
%   menu item and toolbar button click.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1.2.1 $    $Date: 2011/06/07 22:42:56 $


if ishghandle(src, 'uimenu')
    if strcmpi(get(src, 'Checked'), 'on')
        newState = 'off';
    else
        newState = inout;
    end
elseif ishghandle(src, 'uitoggletool')
    newState = get(src, 'State');
    if ~strcmp(newState, 'off')
        newState = inout;
    end
else
    error(message('curvefit:sftoolgui:FitFigure:UnexpectedSourceType'));
end
sftoolgui.sfzoom3d(fitFigure, newState);

end
