function savelineproperties(ds)
%SAVELINEPROPERTIES Save current line properties for later recovery

% $Revision: 1.2.2.2 $  $Date: 2005/03/07 17:24:48 $
% Copyright 2001-2005 The MathWorks, Inc.

lineproperties = {'Color' 'Marker' 'LineStyle' 'LineWidth'};
cml = cell(4,1);
if ~isempty(ds.line) && ishandle(ds.line)
   cml = get(ds.line, lineproperties);
end
ds.ColorMarkerLine = cml;
