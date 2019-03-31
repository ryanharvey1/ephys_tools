%TableOfFitsClient Panel for displaying fits in SFTOOL
%
%   TABLEOFFITSCLIENT

%   Copyright 2008-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $    $Date: 2010/10/25 15:54:50 $

classdef TableOfFitsClient < sftoolgui.Client

    properties (SetAccess = private)
        HFitsManager
        LocalJavaPanel
        FitdevListeners = {};
    end
    
    methods
        function obj = TableOfFitsClient(fitsManager)   
            
            % Get the singleton instance of the client
            obj.JavaClient = javaMethodEDT( 'getInstance', ...
                ['com.mathworks.toolbox.curvefit.' ...
                 'surfacefitting.SFTableOfFitsClient']);
            obj.JavaPanel = javaMethodEDT('getPanel', obj.JavaClient);
            obj.LocalJavaPanel = javaMethodEDT('getPanel', obj.JavaClient);
            
            obj.Name = javaMethodEDT('getName', obj.JavaClient);

            % Same action occurs when fits are "added" or "loaded".
            obj.HFitsManager = fitsManager;
            obj.HFitsManager.addlistener('FitAdded', ...
                @(s, e) obj.fitAdded( e ) );
            obj.HFitsManager.addlistener('FitDeleted', ...
                @(s, e) obj.fitDeleted( e ) );
            obj.HFitsManager.addlistener('FitLoaded', ...
                @(s, e) obj.fitAdded( e ) );
             
            % Setup listeners on this panel
            obj.InternalListeners = [
                addlistener( obj.JavaPanel, 'duplicateFit', ...
                    @(s, e) obj.duplicateFitAction( e ) );
                addlistener( obj.JavaPanel, 'deleteFit', ...
                    @(s, e) obj.deleteFitAction( e ) );
                addlistener( obj.JavaPanel, 'saveToWorkSpace', ...
                    @(s, e) obj.saveToWorkspaceAction( e ) );
                addlistener( obj.JavaPanel, 'selectFit', ...
                    @(s, e) obj.selectFitAction( e ) );
            ];           
        end
             
        function editFitAction(obj, evt)
            name = char( evt.getFitName );
            obj.HFitsManager.editFit(name );
        end
        
        function duplicateFitAction(obj, evt)
            duplicateFit(obj.HFitsManager, evt.getFitUUID );
        end

        function deleteFitAction(obj, evt)
            deleteFit(obj.HFitsManager, evt.getFitUUID );
        end
    
        function saveToWorkspaceAction(obj, evt)
            obj.HFitsManager.saveToWorkspace(evt.getFitUUID);
        end
        
        function selectFitAction(obj, evt)
            selectFit(obj.HFitsManager, evt.getFitUUID );
        end
      
        function fittingDataUpdated(obj, evt)
            javaMethodEDT('clearFittingData', obj.JavaPanel, ...
                evt.HFitdev.FitName, ...
                evt.HFitdev.FittingData.Name, ...
                evt.HFitdev.FitTypeString, ...
                evt.HFitdev.FitState, ...
                evt.HFitdev.FitID);
            javaMethodEDT('clearValidationData', obj.JavaPanel, ...
                evt.HFitdev.FitID);
        end
        
        function validationDataUpdated(obj, evt)
            javaMethodEDT('updateValidationName', obj.JavaPanel, ...
                evt.HFitdev.ValidationData.Name, evt.HFitdev.FitID);
            javaMethodEDT('clearValidationData', obj.JavaPanel, ...
                evt.HFitdev.FitID);
            iUpdateTable(obj, evt.HFitdev);
        end
        
        function fitTypeFitValuesUpdated(obj, evt)
            iUpdateTable(obj, evt.HFitdev);
        end
        
        function fitNameUpdated(obj, evt)
            javaMethodEDT('updateFitName', obj.JavaPanel, ...
                evt.HFitdev.FitName, evt.HFitdev.FitID);
        end
          
        function fitFitted( obj, evt )
            iUpdateTable(obj, evt.HFitdev);
        end

        function obj = fitAdded( obj, evt )
            obj = iAddFitToTable(obj, evt.HFitdev);
            iUpdateTable(obj, evt.HFitdev);
        end     
       
        function obj = fitDeleted( obj, evt )
            javaMethodEDT('removeFitFromTable', obj.JavaPanel, ...
                evt.HFitdev.FitID );
            updateFitdevListeners( obj );
        end
        
        function tableConfig = saveSession(obj)
            tableConfig = sftoolgui.TableOfFitsConfiguration();
            % Get the desktop object
            dt = javaMethodEDT( 'getDesktop', ...
                'com.mathworks.mlservices.MatlabDesktopServices' );
            % Get the Table of Fits desktop client
            dtTOFClient = dt.getClient(obj.Name);
            tableConfig.Visible = dt.isClientShowing(obj.Name);
            % TODO we may want to save location as a preference instead of
            % a configuration element 
            clientLocation = dt.getClientLocation(dtTOFClient);
            %TODO if table is not visible, then path is empty; need a 
            % different way to get that information. For now we won't set
            % the location if it is not empty.
            path = get(clientLocation, 'Path');
            tableConfig.Location = path;
        end
  
         function loadSession(obj, tableConfig)
            % Get the desktop object
            dt = javaMethodEDT( 'getDesktop', ...
                'com.mathworks.mlservices.MatlabDesktopServices' );
            % TODO we may want to save location as a preference instead of
            % a configuration element 
            % Create a desktop location object
            if ~isempty(tableConfig.Location)
                dtLocation = javaMethodEDT( 'create', ...
                    'com.mathworks.widgets.desk.DTLocation', ...
                    tableConfig.Location);
                dt.setClientLocation(obj.Name, dtLocation);
            end
            if tableConfig.Visible
                dt.showClient(obj.Name);
            else
                dt.hideClient(obj.Name);
            end
         end
        
         function closeClient(obj)
             javaMethodEDT('close', obj.JavaClient);
         end
         
         function cleanup(obj)
            javaMethodEDT('clearTable', obj.LocalJavaPanel);
            javaMethodEDT('cleanup',  obj.LocalJavaPanel);
            closeClient(obj);  
         end
    end
    
    methods(Access = private)
        function updateFitdevListeners( obj )
            % Change the source of the Fitdev Listeners to match the list of Fits in the
            % FitsManager. If there are no listeners, then create new ones. If there are no
            % fits, then delete all the listeners.
            
            hFitdevs = obj.HFitsManager.Fits;
            nFits = length( hFitdevs );
            
            if nFits == 0
                % delete all listeners
                obj.FitdevListeners = {};
            else
                % Need to recreate listeners
                obj.FitdevListeners = iCreateFitdevListeners( obj, hFitdevs );                
            end
        end
    end
end

function obj = iAddFitToTable(obj, hFitdev)
   javaMethodEDT('addFit', obj.JavaPanel, ...
                           hFitdev.FitName, ...
                           hFitdev.FitID, ...
                           hFitdev.FitState, ...
                           hFitdev.FitTypeString);
   updateFitdevListeners( obj );
end

function iUpdateTable(obj, HFitdev)

    if HFitdev.FitState == com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.GOOD || ...
       HFitdev.FitState == com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.WARNING  % good or warning fit
        javaMethodEDT('updateFit', obj.JavaPanel, ...
            HFitdev.FitName, ...
            HFitdev.FittingData.Name, ...
            HFitdev.FitTypeString, ...
            HFitdev.Goodness.sse, ...
            HFitdev.Goodness.rsquare, ...
            HFitdev.Goodness.dfe, ...
            HFitdev.Goodness.adjrsquare, ...
            HFitdev.Goodness.rmse, ...
            HFitdev.Output.numparam, ...
            HFitdev.FitID, ...
            HFitdev.FitState);
        javaMethodEDT('updateValidationName', obj.JavaPanel, ...
                    HFitdev.ValidationData.Name, HFitdev.FitID);
        if ~isempty(HFitdev.ValidationGoodness.sse) && ...
           ~isempty(HFitdev.ValidationGoodness.rmse)
               javaMethodEDT('updateValidationData', obj.JavaPanel, ...
                   HFitdev.ValidationGoodness.sse, ...
                   HFitdev.ValidationGoodness.rmse, ...
                   HFitdev.FitID);
        else        
               javaMethodEDT('clearValidationData', obj.JavaPanel, ...
                   HFitdev.FitID);
        end
    else % bad or incomplete
        javaMethodEDT('clearFittingData', obj.JavaPanel, ...
            HFitdev.FitName, ...
            HFitdev.FittingData.Name, ...
            HFitdev.FitTypeString, ...
            HFitdev.FitState, ...
            HFitdev.FitID);
        javaMethodEDT('updateValidationName', obj.JavaPanel, ...
                    HFitdev.ValidationData.Name, HFitdev.FitID);
        javaMethodEDT('clearValidationData', obj.JavaPanel, HFitdev.FitID);
    end
end

function listeners = iCreateFitdevListeners( obj, hFitdevs )
listeners = {
    event.listener( hFitdevs, 'FitNameUpdated',          @(s, e) obj.fitNameUpdated( e ));
    event.listener( hFitdevs, 'FitFitted',               @(s, e) obj.fitFitted( e ));
    event.listener( hFitdevs, 'FittingDataUpdated',      @(s, e) obj.fittingDataUpdated( e ));
    event.listener( hFitdevs, 'ValidationDataUpdated',   @(s, e) obj.validationDataUpdated( e ));
    event.listener( hFitdevs, 'FitTypeFitValuesUpdated', @(s, e) obj.fitTypeFitValuesUpdated( e ));
    };
end

