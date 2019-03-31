function newState = getMenuToggleToolState(obj)
%getMenuToggleToolState gets the new state of a uimenu or uitoggletool

% getMenuToggleToolState(OBJ) gets the state of OBJ which is expected to be
% either a uimenu or a uitoggletool.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/03/22 18:30:10 $

if ishghandle(obj, 'uimenu')
    if strcmpi(get(obj, 'Checked'), 'on')
        newState = 'off';
    else
        newState = 'on';
    end
elseif ishghandle(obj, 'uitoggletool')
    newState = get(obj, 'State');
else
    error(message('curvefit:sftoolgui:getMenuToggleToolState:UnexpectedObjectType'));
end
end