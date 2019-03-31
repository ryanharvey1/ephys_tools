function cfadjustlayout(cffig,showctrl)
%ADJUSTLAYOUT Adjust layout of buttons and graph in figure window

%   $Revision: 1.13.2.5 $  $Date: 2008/10/31 05:55:59 $
%   Copyright 2001-2008 The MathWorks, Inc.

% Get some measurements
ufig = get(cffig,'Units');
fpos = hgconvertunits(cffig,get(cffig,'position'),ufig,'points',0);
fwidth = fpos(3);
fheight = max(1,fpos(4));

% Adjust all button positions
% *** "tags" also defined in cfaddbuttons and cfgetfigurepos, and they must match
tags = {'cfimport' 'cffit' 'cfexclude' 'cfplot' 'cfanalyze'};
hbuttons = findobj(cffig,'Type','uicontrol','Style','pushbutton');
nbuttons = length(tags);
if isempty(hbuttons)
   bheight = 0;
   border = 0;
else
   hb = zeros(nbuttons,1);
   set(hbuttons,'Units','points');
   extents = zeros(nbuttons,4);
   for j=1:nbuttons
      hb(j) = findobj(hbuttons,'flat','Tag',tags{j});
      extents(j,:) = get(hb(j),'Extent');
   end
   bheight = 1.5 * extents(1,4);  % 1.5 * text height
   border = bheight + bheight/2;   % border above buttons
   gutter = bheight/4;            % between buttons
   margin = bheight/2;            % around text within button
   bwidth = extents(:,3)' + 2*margin;
   totalwidth = sum(bwidth) + gutter*(nbuttons-1);
   startpos = max(0, (fwidth/2) - (totalwidth/2));
   bleft = startpos + [0, cumsum(bwidth + gutter)];
   for j=1:nbuttons
      pos = [bleft(j), max(1,fheight-bheight-border), bwidth(j), bheight];
      set(hb(j),'Units','points', 'Position',pos);
   end
end

% Position the axes in the remaining area
ax1 = findall(cffig,'Type','axes','Tag','main');
ax2 = findall(cffig,'Type','axes','Tag','resid');
ytop = max(2,fheight-bheight-border)/fheight;
if isempty(ax2)
   set(get(ax1,'Parent'),'Position',[0 0 1 ytop]); % position uicontainer
else
   set(get(ax1,'Parent'),'Position',[0, ytop/2, 1, ytop/2]);
   set(get(ax2,'Parent'),'Position',[0, 0, 1, ytop/2]);
end

if nargin<2
   showctrl = cfgetset('showaxlimctrl');
end
if isequal(showctrl,'on')
   cfaxlimctrl(cffig,'off');
   cfaxlimctrl(cffig,'on');
end
