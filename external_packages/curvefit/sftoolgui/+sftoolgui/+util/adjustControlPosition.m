function pos = adjustControlPosition(hControl,pos)
%adjustControlPosition Calculate a position that puts a control in the right place
%
%   adjustControlPosition(hControl, pos) adjusts the given position
%   rectangle so that when applied to hControl, hControl will end up in the
%   right place.

%    Copyright 2011 The MathWorks, Inc.
%    $Revision: 1.1.6.1 $    $Date: 2011/05/09 00:39:59 $

if graphicsversion(hControl, 'handlegraphics')
    hFig = ancestor(hControl, 'figure');
    rend = get(hFig, 'Renderer');
    isOpenGL = strcmpi(rend, 'opengl');
    
    if isOpenGL && ishghandle(hControl, 'uipanel')
        % Panels are offset by 1 pixel downwards in OpenGL so add a pixel
        % correction
        pos(2) = pos(2) + 1;
    elseif ishghandle(hControl, 'hgjavacomponent')
        % Apply accumulation of all uipanel parents' border widths
        allW = iGetAllParentBorderWidths(hControl);
        pos(1:2) = pos(1:2) + allW;
    else
        % No correction required for axes, uicontrols
    end
    
    if isOpenGL && ishghandle(get(hControl, 'Parent'), 'uipanel')
        % The offset we added for a parent panel needs to be undone for
        % its contents
        pos(2) = pos(2) - 1;
    end
    
else
    if ~ishghandle(hControl, 'axes')
        % Objects in panels need to be shifted upwards.
        hParent = get(hControl, 'Parent');
        if ishghandle(hParent, 'uipanel')
            W = sftoolgui.util.getPanelBorderWidth(hParent);
            pos(2) = pos(2) + 2*W;
        end
    end
end



function W = iGetAllParentBorderWidths(hControl)

W = 0;
hParent = get(hControl, 'Parent');
while ~isempty(hParent)
    if ishghandle(hParent, 'uipanel')
        W = W + sftoolgui.util.getPanelBorderWidth(hParent);
    end
    hParent = get(hParent, 'Parent');
end
