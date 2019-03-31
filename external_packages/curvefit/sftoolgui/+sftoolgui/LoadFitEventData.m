classdef LoadFitEventData < event.EventData
    %LOADFITEVENTDATA Data for a load fit event
    %
    %   SFTOOLGUI.LOADFITTEVENTDATA
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2008/10/08 17:04:24 $

    properties(SetAccess = 'private', GetAccess = 'public')
         HFitdev;
         HFitFigureConfig;
    end

    methods
        function this = LoadFitEventData( fdev, config )
            this.HFitdev = fdev;
            this.HFitFigureConfig = config;
        end
    end
end
