classdef CopyFitEventData < sftoolgui.FitEventData
    %COPYFITEVENTDATA Data for a FitsManager copy fit event
    %
    %   SURFTOOLGUI.COPYFITTEVENTDATA
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2008/09/13 06:49:32 $

    properties(SetAccess = 'private', GetAccess = 'public')
         HSrcFitdev;
    end

    methods
          function this = CopyFitEventData( fitdev, srcFitdev )
                this = this@sftoolgui.FitEventData(fitdev);
                this.HSrcFitdev = srcFitdev;
          end
     end
end
