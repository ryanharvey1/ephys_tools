function fpos = cfgetfigurepos(cffig,units)
%CFGETFIGUREPOS Get position for a curve fitting figure without uicontrol buttons

%   $Revision: 1.1.6.1 $  $Date: 2004/07/05 17:00:48 $
%   Copyright 2001-2004 The MathWorks, Inc.

% Get some measurements
oldu = get(cffig,'Units');
fpos = get(cffig,'Position');
fpos = hgconvertunits(cffig,fpos,oldu,units,cffig);

% Adjust all button positions
% *** "tags" also defined in cfadjustlayout and cfaddbuttons, and they must match
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
      extents(j,:) = hgconvertunits(cffig,get(hb(j),'Extent'),...
                                    get(hb(j),'Units'),units,cffig);
   end
   bheight = 1.5 * extents(1,4);  % 1.5 * text height
   border = bheight/2;            % border above buttons
end

fpos(4) = fpos(4) - bheight - 2*border;
