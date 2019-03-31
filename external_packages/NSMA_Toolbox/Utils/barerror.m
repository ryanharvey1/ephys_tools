function hh = barerror(Xvals,Yvals,Yerr,XLabels,legendstrings)
%
% barerror(Xvals,Yvals,Yerr [,XLabels,legendstrings])
% 
% plot a grouped bar graph with errorbars
% PL 2000

hh = bar(Xvals,Yvals);
ng = length(hh);                % number of groups of bars

%fh = figure;
%fh2 = bar(Xvals,Yvals, 0.8);
if exist('legendstrings','var')
   legend(legendstrings);
end
hold on;
for ig = 1:ng
   %xc = mean(get(hh(ig),'XData'));  % prior to matlab 7
   xc = mean(get(get(hh(ig),'Children'),'XData'));  % extract XData of the patches of each barseries and get mean position of each patch
   errorbar(xc,Yvals(:,ig),Yerr(:,ig),'dg');
end%for

if exist('XLabels','var')
   set(gca,'XTickLabel',XLabels);
end%if   

hold off;
