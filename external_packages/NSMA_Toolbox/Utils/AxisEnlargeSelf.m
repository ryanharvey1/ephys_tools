function AxisEnlargeSelf

% AxisEnlargeSelf  Allows a subplot to be clicked on to be expanded
%
% AxisEnlargeSelf
%
% INPUTS:
%       (none)
% OUTPUTS:
%       (none)
%
% If you set the 'ButtonDownFcn' property of a subplot to this function,
% then clicking on the subplot will open a new figure with this subplot
% filling the entire figure.  e.g.,
%    h = subplot(nL, ppL, iSP);
%    set(h, 'ButtonDownFcn', 'AxisEnlargeSelf');
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status PROMOTED


axesID = gca;
data = get(gca, 'Children');
titleString = get(get(gca, 'Title'), 'String');
xlabelString = get(get(gca, 'Xlabel'), 'String');
ylabelString = get(get(gca, 'Ylabel'), 'String');

figID = figure;
for iD = length(data):-1:1
  d0 = data(iD);
  if (strcmp(get(d0, 'Type'),'line') == 1)
    x = get(d0, 'XData');
    y = get(d0, 'YData');
    hold on
    plot(x,y, ...
	'color', get(d0, 'Color'), ...
	'marker', get(d0, 'Marker'), ...
	'linestyle', get(d0, 'LineStyle'));
    hold off
  end
end
title(titleString);
xlabel(xlabelString);
ylabel(ylabelString);

zoom on