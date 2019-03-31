classdef FittingConfiguration < sftoolgui.Configuration
    %FITTINGCONFIGURATION configuration file for sftool Fitting panel 
    %
    %   SFTOOLGUI.FITTINGCONFIGURATION
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2008/09/13 06:49:42 $

    properties(SetAccess = 'private', GetAccess = 'private')
         Version = 1;
    end
    
    properties(SetAccess = 'public', GetAccess = 'public')
         Visible = 'on';
    end

    methods
          function this = FittingConfiguration()
          end
    end
end
