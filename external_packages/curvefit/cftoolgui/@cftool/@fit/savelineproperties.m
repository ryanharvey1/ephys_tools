function savelineproperties(fit)
%SAVELINEPROPERTIES Save current line properties for later recovery

% $Revision: 1.2.2.2 $  $Date: 2005/03/07 17:24:53 $
% Copyright 2001-2005 The MathWorks, Inc.

cml = cell(8,1);
lineproperties = {'Color' 'Marker' 'LineStyle' 'LineWidth'};

if ~isempty(fit.line) && ishandle(fit.line)
   cml(1:4) = get(fit.line, lineproperties);
end
if ~isempty(fit.rline) && ishandle(fit.rline)
   cml(5:8) = get(fit.rline, lineproperties);
end

fit.ColorMarkerLine = cml;