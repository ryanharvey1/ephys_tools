classdef HideablePanel <sftoolgui.Panel
    %HIDEABLEPANEL A panel that can be used by An HG figure for use with SFTOOL
    %
    %   HIDEABLEPANEL(fig)

    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2008/09/13 06:49:44 $ 

    events
        PanelClosing;
    end
    
    % TODO add a small "<< or >>" in the top right corner to be used to hide the
    % panel
    methods
        function this = HideablePanel(parentFig)
            this = this@sftoolgui.Panel(parentFig);
        end
     end
 end


