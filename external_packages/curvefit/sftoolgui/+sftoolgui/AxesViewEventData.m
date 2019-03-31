classdef AxesViewEventData < event.EventData
    %AXESVIEWEVENTDATA Data for a AxesViewModel event
    %
    %   SFTOOLGUI.AXESVIEWEVENTDATA
    
    %   Copyright 2011 The MathWorks, Inc.
    %   $Revision: 1.1.4.1 $    $Date: 2011/07/18 00:30:47 $
    
    properties(SetAccess = 'private', GetAccess = 'public')
        AxesViewModel;
    end
    
    methods
        function this = AxesViewEventData( axesViewModel )
            this.AxesViewModel = axesViewModel;
        end
    end
end
