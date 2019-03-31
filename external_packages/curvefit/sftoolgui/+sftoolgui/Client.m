%CLIENT Abstract client for use in Surface Fitting Tool.
%
%   CLIENT

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $    $Date: 2009/03/20 18:17:12 $ 

classdef Client < handle
    properties(SetAccess = 'protected', GetAccess = 'protected')
        % JavaPanel --  The Java version of this panel. 
        JavaPanel;
        JavaClient;
        % InternalListeners -- Listeners to java events that are internal to the
        % panel
        InternalListeners;
    end

    properties(SetAccess = 'protected', GetAccess = 'public')
        % Name -- The name of the panel - this should not be translated
        Name;
    end
    
    methods
        function obj = Client
        end
        
        function addClient( obj, dt, dtl )
            % ADDCLIENT -- Add the panel to the desktop
            %     ADDCLIENT( HPANEL, DT, DTL ) adds the panel HPANEL to the
            %     MATLAB desktop DT in the position DTL (DTLocation)
            %
            % TODO: Is this the right place to be positioning the panel?
            % Then we position the client in the desktop
            %
            % 'Name' is actually controlled by setClientName and setTitle
            % in the JavaClient
            javaMethodEDT( 'addClient', dt, obj.JavaClient, obj.Name, true, dtl, true );
        end
    end

end
