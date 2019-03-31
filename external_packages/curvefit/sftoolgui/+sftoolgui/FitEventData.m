classdef FitEventData < event.EventData
    %FITEVENTDATA Data for a FitsManager event
    %
    %   SURFTOOLGUI.FITTEVENTDATA
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2008/09/13 06:49:37 $

    properties(SetAccess = 'private', GetAccess = 'public')
           HFitdev;
    end

    methods
        function this = FitEventData( fitdev )
                this.HFitdev = fitdev;
        end
    end
end
