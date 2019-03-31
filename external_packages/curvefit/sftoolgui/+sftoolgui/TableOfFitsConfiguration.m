classdef TableOfFitsConfiguration < sftoolgui.Configuration
    % TABLEOFFITSCONFIGURATON configuration file for sftool TableOfFits 
    
    %   SFTOOLGUI.TABLEOFFITSCONFIGURATION
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.2 $    $Date: 2008/10/31 05:58:01 $

    properties(SetAccess = 'private', GetAccess = 'private')
         Version = 1;
    end
    properties(SetAccess = 'public', GetAccess = 'public')
         Visible = true;
         Location = 'S';
    end

    methods
          function this = TableOfFitsConfiguration()
          end
     end
end
