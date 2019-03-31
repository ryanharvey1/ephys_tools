classdef (Sealed) AxesViewModel < handle
    %AxesViewModel is the model for axes View and Limits properties
    
    %   Copyright 2011 The MathWorks, Inc.
    %   $Revision: 1.1.4.2 $    $Date: 2011/08/01 00:51:23 $
    
    events (NotifyAccess = private)
        %ThreeDViewAngleChanged-- fired when the ThreeDViewAngle property
        %changes
        ThreeDViewAngleChanged
        %LimitsChanged-- fired when any of the Limit properties change
        LimitsChanged
    end
    
    properties( SetAccess = 'private', GetAccess = 'public')
        
        % XInputLimits are the limits for the x input
        XInputLimits = [-1, 1] ;
        
        % YInputLimits are the limits for the y input
        YInputLimits = [-1, 1] ;
        
        % ResponseLimits are the limits for the response data
        ResponseLimits = [-1, 1] ;
        
        % ResidualLimits are the limits for the residuals response data
        ResidualLimits = [-1 1];
    end
    
    properties( SetAccess = 'public', GetAccess = 'public', Dependent)
        % ThreeDViewAngle is view angle of a non-curve plot
        ThreeDViewAngle = sftoolgui.util.DefaultViewAngle.ThreeD;
    end
    
    properties( SetAccess = 'private', GetAccess = 'private')
        % PrivateThreeDViewAngle is view angle of a non-curve plot that is
        % associated with the public dependent ThreeDViewAngle property
        PrivateThreeDViewAngle = sftoolgui.util.DefaultViewAngle.ThreeD;
    end
    
    methods
        function set.ThreeDViewAngle(this, view)
            this.PrivateThreeDViewAngle = view;
            notify (this, 'ThreeDViewAngleChanged', sftoolgui.AxesViewEventData( this));
        end
        
        function viewAngle = get.ThreeDViewAngle(this)
            viewAngle = this.PrivateThreeDViewAngle;
        end
        
        function setLimits(this, input, output, residuals)
            % setLimits will set all limit properties, so that a
            % LimitsChanged event will include up-to-date values for all
            % the limits.
            setInputs(this, input);
            setResponseLimits(this, output)
            setResidualLimits(this, residuals)
            notify (this, 'LimitsChanged', sftoolgui.AxesViewEventData( this));
        end
        
        function setLimit(this, property, limit)
            % setLimit will set the LIMIT of the specified PROPERTY
            setALimit(this, property, limit);
            notify (this, 'LimitsChanged', sftoolgui.AxesViewEventData( this));
        end
    end
    
    methods(Access = private)
        function setInputs(this, inputs)
            % setInputs sets XInputLimits and YInputLimits.
            setALimit( this, 'XInputLimits', inputs{1} );

            % If there is only one input limit ...
            if length(inputs) == 1  
                % ... then we assume that this is for curve data and we set
                % the YInputLimits to the default value.
                this.YInputLimits = [-1 1];
            else 
                % ... otherwise we set the YInputLimits to the requested value.
                setALimit( this, 'YInputLimits', inputs{2} );
            end
        end
        
        function setResponseLimits(this, responseLimits)
            % setResponseLimits will set the ResponseLimits property if
            % responseLimits is not empty.
            setALimit( this, 'ResponseLimits', responseLimits );
        end
        
        function setResidualLimits(this, residualLimits)
            % setResidualLimits will set the ResidualLimits property if
            % residualLimits is not empty.
            setALimit( this, 'ResidualLimits', residualLimits );
        end
        
        function setALimit( this, name, value )
            % setALimit   Set a single limit.
            %
            % This method does not fire an event and has protection against
            % invalid values for the limit.
            
            % Ignore limits in the following cases:
            %     value is empty
            %     value is not 1x2
            %     value(1) is not less than value(2) 
            %     value is not finite.
            if isempty( value ) || ...
                ~isequal( size( value ), [1, 2] ) || ...
                diff( value ) <= 0 || ...
                ~all( isfinite( value ) );
                % Ignore limits
            else
                % The given value has passed all our tests. We should set
                % the appropriate property.
                this.(name) = value;
            end
        end
    end
end


