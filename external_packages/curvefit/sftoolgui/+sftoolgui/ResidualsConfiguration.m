classdef ResidualsConfiguration < sftoolgui.Configuration
    %RESIDUALSCONFIGURATION configuration file for sftool Residuals plot 
    %
    %   SFTOOLGUI.RESIDUALSCONFIGURATION
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2008/09/13 06:49:47 $

    properties(SetAccess = 'private', GetAccess = 'private')
         Version = 1;
    end

    properties(SetAccess = 'public', GetAccess = 'public')
         Visible = 'off';
    end
    
    methods
          function this = ResidualsConfiguration()
          end
    end
end
