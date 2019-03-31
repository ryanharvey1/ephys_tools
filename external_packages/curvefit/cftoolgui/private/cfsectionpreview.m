function imsource = cfsectionpreview(outlier, width, height)
% CFSECTIONPREVIEW For use by CFTOOL

%   $Revision: 1.2.2.4 $
%   Copyright 2001-2010 The MathWorks, Inc.

if nargin < 2
    width = 180;
    height = 180;
end

imsource = cftoolimagefromplot( @nPlotSection, width, height );

    function hPlot = nPlotSection( ax )
        hPlot = [];
        
        xlim = [0 4];
        ylim = [0 4];
        set( ax, 'Position',[.05 .05 .9 .9], ...
            'XTick',[],'YTick',[], ...
            'Box','on', ...
            'Visible','off', 'XLim',xlim,'YLim',ylim);
        
        gr = [.9 .9 .9];
        o = handle(outlier);
        
        xlo = o.domainLow;
        if ~isempty(xlo)
            hPlot(end+1) = patch([0 1 1 0], [0 0 4 4], gr,'LineStyle','none','Parent',ax);
        end
        
        xhi = o.domainHigh;
        if ~isempty(xhi)
            hPlot(end+1) = patch([3 4 4 3], [0 0 4 4], gr,'LineStyle','none','Parent',ax);
        end
        
        ylo = o.rangeLow;
        if ~isempty(ylo)
            hPlot(end+1) = patch([0 4 4 0], [0 0 1 1], gr,'LineStyle','none','Parent',ax);
        end
        
        yhi = o.rangeHigh;
        if ~isempty(yhi)
            hPlot(end+1) = patch([0 4 4 0], [3 3 4 4], gr,'LineStyle','none','Parent',ax);
        end
    end
end
