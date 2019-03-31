%M Wrapper for SFGroup.java
%
%   GROUP

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $    $Date: 2009/04/15 22:56:36 $

classdef Group < handle

    properties (SetAccess = private)
        SFTool
        HFitsManager
        GroupBase
        FitdevListeners = {};
        InternalListeners ;
    end
    
    methods
        function this = Group(sftool)
           
            this.SFTool = sftool;
            % Get the singleton instance of the group
            this.GroupBase = javaMethodEDT( 'getInstance', ...
                'com.mathworks.toolbox.curvefit.surfacefitting.SFGroup');
            
            this.HFitsManager = sftool.HFitsManager;
            
             % Same action occurs when fits are "added" or "loaded".
            this.HFitsManager.addlistener('FitAdded', ...
                              @(s, e) this.fitAdded( e ) );
            this.HFitsManager.addlistener('FitLoaded', ...
                              @(s, e) this.fitAdded( e ) );
            this.HFitsManager.addlistener('FitDeleted', ...
                              @(s, e) this.fitDeleted( e ) );
             
            % Setup listeners on group events
            this.InternalListeners = [
                  addlistener( this.GroupBase, 'newCrossSectionPlot',...
                     @(s, e) this.SFTool.newCrossSectionPlotAction( e ) );
                  addlistener( this.GroupBase, 'figureActivated', ...
                     @(s, e) this.SFTool.figureActivatedAction( e ) );
                  addlistener( this.GroupBase, 'newFit', ...
                     @(s, e) this.SFTool.newFitAction( ) );
                  addlistener( this.GroupBase, 'closeTool', ...
                     @(s, e) this.SFTool.closeSftool( ) );
                  addlistener( this.GroupBase, 'groupClosed', ...
                     @(s, e) this.SFTool.closeAction( ) );
                  addlistener( this.GroupBase, 'groupClosing', ...
                     @(s, e) this.SFTool.closingAction( ) );
                  addlistener( this.GroupBase, 'selectFit', ...
                     @(s, e) this.selectFitAction( e ) );
                  addlistener( this.GroupBase, 'loadSession', ...
                     @(s, e) this.SFTool.session( 'load' ) );
                  addlistener( this.GroupBase, 'saveSession', ...
                     @(s, e) this.SFTool.session( 'save' ) );
                  addlistener( this.GroupBase, 'saveSessionAs', ...
                     @(s, e) this.SFTool.sessionSaveAs() );
                  addlistener( this.GroupBase, 'clearSession', ...
                     @(s, e) this.SFTool.session( 'clear' ) );
                  addlistener( this.GroupBase, 'sftoolHelp', ...
                     @(s, e) this.SFTool.sftoolHelp );
                  addlistener( this.GroupBase, 'cftoolHelp', ...
                     @(s, e) this.SFTool.cftoolHelp );
                  addlistener( this.GroupBase, 'demosHelp', ...
                     @(s, e) this.SFTool.demosHelp );
                  addlistener( this.GroupBase, 'aboutHelp', ...
                     @(s, e) this.SFTool.aboutHelp );
                  addlistener( this.GroupBase, 'generateMFile', ...
                     @(s, e) this.generateMFileAction())
               ];
        
        end
        
        function addGroup(this, dt)
            dt.addGroup(this.GroupBase);
        end
        
        function setDefaultLocation(this)
            javaMethodEDT( 'setDefaultLocation', this.GroupBase);
        end
        
        function setGroupTitle(this)
            if isempty(this.SFTool.SessionName)
                title = this.SFTool.GroupName;
            else
                [p, name] = fileparts(this.SFTool.SessionName);
                title = sprintf('%s - %s', this.SFTool.GroupName, name);
            end
            javaMethodEDT('setGroupTitle', this.GroupBase, title);
        end
        
          function approveClose(this)
             javaMethodEDT('approveGroupClose', this.GroupBase);
        end
        
        function vetoClose(this)
             javaMethodEDT('vetoGroupClose', this.GroupBase);
        end
        
        function generateMFileAction( this)
            sftoolgui.generateMCode( this.SFTool );
        end
        
        function selectFitAction(this, evt)
            selectFit(this.HFitsManager, evt.getFitUUID );
        end
        
        function fitNameUpdated(this, evt)
            javaMethodEDT('updateFitName', this.GroupBase, ...
                                           evt.HFitdev.FitName, ...
                                           evt.HFitdev.FitID);
        end
         
        function this = fitAdded(this, evt )
            javaMethodEDT('addFit', this.GroupBase, ...
                           evt.HFitdev.FitName, ...
                           evt.HFitdev.FitID);
            updateFitdevListeners( this );
        end
        
        function this = fitDeleted( this, evt )
            javaMethodEDT('deleteFit', this.GroupBase, evt.HFitdev.FitID );
            updateFitdevListeners( this );
        end
        
        function deleteGroup(this)
            delete(this.InternalListeners);
            delete(this)
        end
    end
    
    methods(Access = private )
        function updateFitdevListeners( this )
            % Change the source of the Fitdev Listeners to match the list of
            % Fitdevs in the Fits Manager. If there are no listeners, then
            % create new ones. If there are no  fits, then delete the listeners.
            hFitdevs = this.HFitsManager.Fits;
            
            if isempty( hFitdevs );
                % delete all listeners
                this.FitdevListeners = {};
                
            else
                % Need to create listeners
                this.FitdevListeners = event.listener( hFitdevs, 'FitNameUpdated', ...
                    @(s, e) this.fitNameUpdated(e) );
            end
        end
        
    end
end

