function H = ParentFigureHandle(objHandle)

% ParentFigureHandle  Gets figure handle above current object
%
% H = ParentFigureHandle(objHandle)
%
% INPUTS:
%       objHandle - valid object handle
% OUTPUTS:
%       H - figure handle containing object or root if no figure found
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status PROMOTED


H = objHandle;
while ~strcmp(get(H, 'Type'),'figure') & ~strcmp(get(H, 'Type'),'root')
   H = get(H, 'Parent');
end