classdef Session
    % SESSION for sftool 
    %
    %   SFTOOLGUI.SESSION
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2008/10/31 05:57:57 $

    properties(SetAccess = 'private', GetAccess = 'private')
        Version = 1;
    end    
      
    properties(SetAccess = 'private', GetAccess = 'public')
        CurveFittingToolboxVersion = ver('curvefit');
        AllFitdevsAndConfigs ;
        TableConfig ;
    end

    methods
        function this = Session(allFitdevsAndConfigs, tableConfig)
            this.AllFitdevsAndConfigs = allFitdevsAndConfigs;
            this.TableConfig = tableConfig;
        end
    end
end
