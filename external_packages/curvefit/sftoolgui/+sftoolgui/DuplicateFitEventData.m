classdef DuplicateFitEventData < event.EventData
    %DUPLICATEFITEVENTDATA Data for a FitsManager duplicate fit event
    %
    %   SFTOOLGUI.DUPLICATEFITTEVENTDATA
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.2 $    $Date: 2008/10/31 05:57:42 $

    properties(SetAccess = 'private', GetAccess = 'public')
         SourceFitUUID;
         HFitdev;
    end

    methods
          function this = DuplicateFitEventData( sID,  dF)
                this.SourceFitUUID = sID;
                this.HFitdev = dF; %duplicated fit fitdev
          end
     end
end
