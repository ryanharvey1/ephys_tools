function panModeCallback(fitFigure, src, ~)
% panModeCallback Pan mode callback
%
%   panModeCallback(fitFigure, SRC, EVENT) is the callback to the pan menu
%   item and toolbar button click.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/03/22 18:30:19 $

sftoolgui.sfpan3d(fitFigure, sftoolgui.util.getMenuToggleToolState(src));
end