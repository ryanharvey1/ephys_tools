function hh = barerror2(Xvals,Yvals,Yerr,XLabels)
%
% barerror2(Xvals,Yvals,Yerr [,XLabels])
%
% Xvals ... array of X values
% Yvals ... array of Y values (bar heights)
% Yerr ...  array of errorbar sizes of Y values to be superimposed on bars
% XLabels ... cell array of text strings (same size as X) to label the bars on the x-axis 
% 
% plot a simple (ungrouped) bar graph with errorbars and text labels
% PL 2000

hh = bar(Xvals,Yvals);
ng = length(hh);                % number of groups of bars
if ng ~= 1
    error('Yvals must be an array, not a matrix');
end

%fh = figure;
%fh2 = bar(Xvals,Yvals, 0.8);
if exist('legendstrings','var')
   legend(legendstrings);
end
hold on;
   %xc = mean(get(hh(ig),'XData'));  % prior to matlab 7
xc = mean(get(get(hh(ig),'Children'),'XData'));  % extract XData of the patches of each barseries and get mean position of each patch
errorbar(xc,Yvals,Yerr,Yerr,'dg');

if exist('XLabels','var')
   set(gca,'XTickLabel',XLabels);
end%if   

hold off;
