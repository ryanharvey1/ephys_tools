classdef CloseablePanel <sftoolgui.Panel
    %CLOSEABLEPANEL A panel that can be used with SFTOOL
    %
    %   CLOSEABLEPANEL(parent)

    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.2 $    $Date: 2011/04/11 16:09:31 $      

    events
        PanelClosing;
    end
    
    methods
        function this = CloseablePanel(parent)
            this = this@sftoolgui.Panel(parent);
        end
     end
 end


