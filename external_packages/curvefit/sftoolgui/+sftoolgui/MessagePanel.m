classdef MessagePanel <sftoolgui.Panel
    %MessagePanel A panel that display a message
    %
    %   MessagePanel(parent, message)
    
    %   Copyright 2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.1.2.1 $    $Date: 2011/07/18 00:30:54 $
    
    properties( SetAccess = 'private', GetAccess = 'private')
        % HTextControl is the text uicontrol for the message
        HTextControl
        MessageHeightExtent
    end
    
    properties (Constant, GetAccess = 'private')
        % VerticalFactor is the percentage of the panel's height which is
        % used to set the HTextControl's y position. The position should be
        % a little higher than centered.
        VerticalFactor = .57
    end
    
    methods
        function this = MessagePanel( parent, theMessage )
            % MessagePanel   Construct a MessagePanel.
            %
            %   sftoolgui.MessagePanel( parent, theMessage ) constructs a
            %   MessagePanel. PARENT must be a figure or uipanel.
            %   theMessage is the message that will be displayed in the
            %   panel.
            
            % Make sure parent is either a figure or a uipanel
            if ~(ishghandle(parent, 'figure') || ishghandle(parent, 'uipanel'))
                error(message('curvefit:sftoolgui:MessagePanel:InvalidParent'));
            end
            
            % Make sure theMessage is a string
            if ~ischar(theMessage)
                error(message('curvefit:sftoolgui:MessagePanel:InvalidMessage'));
            end
            
            this = this@sftoolgui.Panel(parent);
            this.Tag = 'sftoolMessageUIPanel';
            this.HTextControl = uicontrol('Parent', this.HUIPanel, ...
                'HorizontalAlignment', 'center', ...
                'Style', 'Text', 'FontSize', 11, 'Units', 'pixels', ...
                'String', theMessage);
            
            extent = get(this.HTextControl, 'Extent');
            this.MessageHeightExtent = extent(4);
            
            this.layoutPanel();
        end
    end
    
    methods(Access = protected)
        function layoutPanel(this)
            % Pass on call to the superclass
            layoutPanel@sftoolgui.Panel(this);
            
            if ~isempty(this.HTextControl)
                innerpos = this.InnerPosition;
                
                idealY = max(1, floor(innerpos(4)*this.VerticalFactor));
                
                % Adding the height extent will insure the text will be
                % shown even when the height of the control gets very small
                controlPosition = [1 1 innerpos(3) idealY + this.MessageHeightExtent];
                
                % Clip the position rectangle to be within innerpos
                controlPosition = sftoolgui.util.clipPosition(innerpos, controlPosition);
                
                % Apply correction
                controlPosition = sftoolgui.util.adjustControlPosition(this.HTextControl, controlPosition);
                
                % Set position
                set(this.HTextControl, 'Position', controlPosition);
            end
        end
    end
end
