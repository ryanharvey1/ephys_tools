classdef PlotLayoutPanel < sftoolgui.Panel
%PLOTLAYOUTPANEL A SFTOOL panel for displaying and laying out plot panels
%
%   PLOTLAYOUTPANEL(hFIG) creates an instance of a PlotLayoutPanel object.
%   This panel class is given handles to a trio of plot panels and
%   positions them correctly.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/05/09 00:39:47 $
    
    % Handles to the three plot panels that need to be positioned
    properties(Access = private)
        SurfacePanel
        ContourPanel
        ResidualsPanel 
    end
    
    properties(Access = private, Transient)      
        % Listener on child visibility that triggers a re-layout
        ChildVisibleListener 
    end
    
    
    methods
        function this = PlotLayoutPanel(parent)
        %PlotLayoutPanel  Constructor for the PlotLayoutPanel class
        %
        %   PlotLayoutPanel(hParent) creates a new instance of the
        %   PlotLayoutPanel class in the specified parent. 
        
            this = this@sftoolgui.Panel(parent);
            
            set(this.HUIPanel, 'BorderType', 'none');
        end
        
        function setPanels(this, hSurface, hResiduals, hContour)
        %setPanels  Set the plot panels to position
        %
        %   setPanels(obj, hSurface, hResiduals, hContour) sets the three
        %   plot panel handles that this object should position.
            
            % Store handles to the individual plot panels
            this.SurfacePanel = hSurface;
            this.ContourPanel = hContour;
            this.ResidualsPanel = hResiduals;
            
            % Add a listener to the Visible property of all of our
            % uipanel's children.  We will redo the layout if any change.
            this.ChildVisibleListener = ...
                curvefit.proplistener(get(this.HUIPanel, 'Children'), 'Visible', @this.updateLayout);
                
            this.layoutPanel();
        end
       
    end
       
    methods(Access = protected)
        function layoutPanel(this)
            %layoutPanel Update the contents of the Panel
            %
            %   layoutPanel is called when the size of the contained uipanel
            %   object changes. This method positions the three plot panels
            %   correctly.
            
            % Pass on call to the superclass
            layoutPanel@sftoolgui.Panel(this);
            
            plotsVisible = struct( ...
                'surfaceVisible', iIsPanelVisible(this.SurfacePanel), ...
                'residualsVisible', iIsPanelVisible(this.ResidualsPanel), ...
                'contourVisible', iIsPanelVisible(this.ContourPanel));
            numVisiblePlots = iNumberOfVisiblePlots(plotsVisible); 
            
            innerpos = this.InnerPosition;
            
            switch numVisiblePlots
                case 1 % one plot
                    positionOnePlot(this, innerpos, plotsVisible);
                case 2 % two plots
                    positionTwoPlots(this, innerpos, plotsVisible);
                case 3 % three plots
                    positionThreePlots(this, innerpos);
                otherwise
                    % Nothing to position
            end
        end
    end
    
    
    
    methods(Access = private)
        function updateLayout(this, ~, ~)
        %updateLayout callback that updates the layout of this panel
            this.layoutPanel();
        end
        
        
        function positionOnePlot(this, innerpos, plotsVisible)
            % positionOnePlot sets the position of one plot if there is only one plot
            % visible.
            hPlotPanel = [];
            if plotsVisible.surfaceVisible
                hPlotPanel = this.SurfacePanel;
            elseif plotsVisible.residualsVisible
                hPlotPanel = this.ResidualsPanel;
            elseif plotsVisible.contourVisible
                hPlotPanel = this.ContourPanel;
            end
            
            % The plot panel is given all of innerpos
            iSetChildPanelPosition(innerpos, hPlotPanel, innerpos);
        end
        
        function positionTwoPlots(this, innerpos, plotsVisible)
            % positionTwoPlots set the position of two plots if exactly two plots are
            % visible.
            
            hTopPanel = [];
            hBottomPanel = [];
            if iAreOnlySurfaceResidualsVisible(plotsVisible)
                hTopPanel = this.SurfacePanel;
                hBottomPanel = this.ResidualsPanel;
            elseif iAreOnlySurfaceContourVisible(plotsVisible)
                hTopPanel = this.SurfacePanel;
                hBottomPanel = this.ContourPanel;
            elseif iAreOnlyResidualsContourVisible(plotsVisible)
                hTopPanel = this.ResidualsPanel;
                hBottomPanel = this.ContourPanel;
            end
            
            bottomH = floor(innerpos(4)/2);
            
            toppos = [innerpos(1), innerpos(2) + bottomH, innerpos(3), innerpos(4) - bottomH];
            bottompos = [innerpos(1), innerpos(2), innerpos(3), bottomH];
            
            iSetChildPanelPosition(innerpos, hTopPanel, toppos);
            iSetChildPanelPosition(innerpos, hBottomPanel, bottompos);
        end
        
        function positionThreePlots(this, innerpos)
            % positionThreePlots sets the positions of three plots if all three plots
            % are visible.
            
            hLeftPanel = this.ContourPanel;
            hTopRightPanel = this.SurfacePanel;
            hBottomRightPanel = this.ResidualsPanel;
            
            leftW = floor(innerpos(3)/2);
            bottomH = floor(innerpos(4)/2);
       
            leftpos = [innerpos(1), innerpos(2), leftW, innerpos(4)];   
            bottomrightpos = [innerpos(1) + leftW, innerpos(2), innerpos(3) - leftW, bottomH];
            toprightpos = [innerpos(1) + leftW, innerpos(2) + bottomH, innerpos(3) - leftW, innerpos(4) - bottomH];
            
            iSetChildPanelPosition(innerpos, hLeftPanel, leftpos);
            iSetChildPanelPosition(innerpos, hTopRightPanel, toprightpos);
            iSetChildPanelPosition(innerpos, hBottomRightPanel, bottomrightpos);
        end

    end
end


function vis = iIsPanelVisible(hPanel)
%iIsPanelVisible Test whether a panel exists and is visible

vis = false;
if ~isempty(hPanel)
    vis = strcmpi(hPanel.Visible, 'on');
end
end


function iSetChildPanelPosition(innerpos, hPanel, pos)
%iSetChildPanelPosition Set a plot panel to a given position

% Clip the position rectangle to be within innerpos
hPanel.Position = sftoolgui.util.clipPosition(innerpos, pos);
end


function num = iNumberOfVisiblePlots(plotsVisible)
%iNumberOfVisiblePlots Return the number of visible plots

num = 0;
if plotsVisible.surfaceVisible
    num = num + 1;
end
if plotsVisible.residualsVisible
    num = num + 1;
end
if plotsVisible.contourVisible
    num = num + 1;
end
end

function tf = iAreOnlySurfaceResidualsVisible(plotsVisible)
% iAreOnlySurfaceResidualsVisible returns true if the surface and residual
% plots are visible and the contour plot is not.
tf = plotsVisible.surfaceVisible &&  plotsVisible.residualsVisible && ~plotsVisible.contourVisible;
end

function tf = iAreOnlySurfaceContourVisible(plotsVisible)
% iAreOnlySurfaceContourVisible returns true if the surface and contour
% plots are visible and the residuals plot is not.
tf = plotsVisible.surfaceVisible && ~plotsVisible.residualsVisible && plotsVisible.contourVisible;
end

function tf = iAreOnlyResidualsContourVisible(plotsVisible)
% iAreOnlyResidualsContourVisible returns true if the residuals and contour
% plots are visible and the surface plot is not
tf = ~plotsVisible.surfaceVisible && plotsVisible.residualsVisible && plotsVisible.contourVisible;
end
