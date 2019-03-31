classdef CloseFitEventData < event.EventData
    %CLOSEFITEVENTDATA Data for a close fit event
    %
    %   SFTOOLGUI.CLOSEFITTEVENTDATA
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2008/09/13 06:49:29 $

    properties(SetAccess = 'private', GetAccess = 'public')
         FitUUID;
         FitFigureConfig;
    end

    methods
          function this = CloseFitEventData( fID, config )
                this.FitUUID = fID;
                this.FitFigureConfig = config;
          end
     end
end
