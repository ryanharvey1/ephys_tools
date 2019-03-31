classdef FitDefinition < handle
    %Surface Fitting FitDefinition holds Fitdev information that differs
    %depending on whether the fit is for curves or surfaces
    
    %   Copyright 2010 The MathWorks, Inc.
    %   $Revision: 1.1.6.2 $    $Date: 2011/05/09 00:39:40 $
    
    properties(SetAccess = 'private', GetAccess = 'private')
        % Version  - class version number
        Version = 1;
    end
    
    properties
        % Fit Type Object
        Type ; 
        % Fit Options
        Options ;
    end  
    
    methods
        function this = FitDefinition(type, options)
            % Construct an sftoolgui.FitDefinition object using type and options if supplied.
            %
            % FitDefinition()
            % Creates a default FitDefinition object
            % FitDefinition(type, options)
            % Creates a FitDefinition using type and options if supplied.
            
            if (nargin == 2)
                if ~isa(type, 'fittype')
                    error(message('curvefit:sftoolgui:FitDefinition:InvalidType'));
                end
                
                if ~isa(options, 'curvefit.basefitoptions')
                    error(message('curvefit:sftoolgui:FitDefinition:InvalidOptions'));
                end
         
                this.Type = type;
                this.Options = options;
            else
                this.Type = fittype( 'linearinterp', 'numindep', 2 );
            	this.Options = fitoptions(this.Type);
                this.Options.Normalize = 'on';
            end          
      
        end
    end
end
        
