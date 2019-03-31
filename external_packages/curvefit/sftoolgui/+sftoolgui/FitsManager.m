classdef FitsManager < handle
    %FITSMANAGER Fits manager for SFTOOL
    %
    %   The FITSMANAGER provides access to fit information
        
    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.8 $    $Date: 2011/05/09 00:39:42 $

    % FitAdded corresponds to new fit buttons on the data panel, edit fit
    % panel and the table of fits panel
    events
        FitCanceled
        FitAdded
        FitDuplicated
        FitDeleted
        FitSavedToWorkspace
        FitSelected
        FitClosed
        FitLoaded
        FitsManagerSessionChanged
        SessionCleared
    end

    properties (GetAccess = public, SetAccess = private)
        Fits = {};
        %TODO where should be keep/control the fit name counter??
        FitNameCount = 1;
        FitdevListeners = {};
    end

    methods
        function this = FitsManager

        end
        
        % The following methods do the actual actions
        % The following are the fit methods
        
        function hFitdev = newFit(this, sftoolguiData)
            % This call seems particularly slow. Is there a m UUID?
            % TODO do we need fitUUIDs with the newest scheme?
            fitUUID = javaMethodEDT('randomUUID', 'java.util.UUID');
            [fitName, this.FitNameCount] = iGetUniqueFitName(this);
            hFitdev = sftoolgui.Fitdev(fitName, fitUUID, sftoolguiData);
            this = addFit(this, hFitdev);
            notify (this, 'FitAdded', sftoolgui.FitEventData( hFitdev ));
        end
        
        % The Fitdev method checks for non empty strings.
        % We want to make an additional check here to ensure names are
        % unique.
        function theMessage = checkFitNameFcn(this, fitName, targetFitdev)
            if isFitNameTaken(this, fitName, targetFitdev)
				theMessage = message('curvefit:sftoolgui:FitsManager:duplicateFitName');
            else
                theMessage = [];
            end
        end
                
        function loadSession(this, allFitdevsAndConfigs)
            for i=1:length(allFitdevsAndConfigs);
                loadFit(this, ...
                    allFitdevsAndConfigs{i}.Fitdev, ...
                    allFitdevsAndConfigs{i}.Config);
            end
        end
        
        function hFitdev = loadFit(this, hFitdev, fitFigureConfig)
            this = addFit(this, hFitdev);
            notify (this, 'FitLoaded', sftoolgui.LoadFitEventData( hFitdev, fitFigureConfig));
        end
               
        function closeFit(this, fitUUID, fitFigureConfig)
            notify (this, 'FitClosed', sftoolgui.CloseFitEventData(fitUUID, fitFigureConfig));
        end
        
        function duplicatedFit = duplicateFit(this, fitUUID)
             hSrcFitdev = getFitFromUUID(this.Fits, fitUUID);
             newName = iGetNameForDuplicatedFit(hSrcFitdev.FitName, this.Fits);
             duplicatedFit = createADuplicate(hSrcFitdev, newName);
             this = addFit(this, duplicatedFit);
             notify (this, 'FitAdded', sftoolgui.FitEventData( duplicatedFit ));
             notify (this, 'FitDuplicated', sftoolgui.DuplicateFitEventData(fitUUID, duplicatedFit)); 
        end
             
        function selectFit(this, fitUUID)
            hFitdev = getFitFromUUID(this.Fits, fitUUID);
            notify (this, 'FitSelected', sftoolgui.FitEventData( hFitdev)); 
        end
        
        function editFit(this, name)
            hFitdev = getFitFromName(this.Fits, name);
            notify (this, 'FitEdited', sftoolgui.FitEventData( hFitdev)); 
        end
        
        function clearSession(this)
            deleteAllFits(this);
            this.FitNameCount = 1;
            notify (this, 'SessionCleared'); 
        end
        
        function deleteAllFits(this)
            n = size(this.Fits, 2);
            for i = 1:n
                 notify (this, 'FitDeleted', sftoolgui.FitEventData( this.Fits{i} )); 
            end
            this.Fits = {};
            this.FitdevListeners = {};
        end
        
        function deleteFit(this, fitUUID)
            [hFitdev, i] = getFitFromUUID(this.Fits, fitUUID);
            notify (this, 'FitDeleted', sftoolgui.FitEventData( hFitdev )); 
            updateFitdevListeners( this );
            this.Fits(i) = [];
        end
        
        function saveToWorkspace(this, fitUUID)
            hFitdev = getFitFromUUID(this.Fits, fitUUID);
            saveFitToWorkspace(hFitdev);
        end
        
        % isFitNameTaken is used in at least two places: when new fits are
        % created and when users change existing fit names. When checking
        % in the case of existing fits, we want to skip that fit when the
        % names are being checked. 
        function taken = isFitNameTaken(this, fitName, targetFitdev)
            n = size(this.Fits, 2);
            taken = false;
            for i = 1:n
                fit = this.Fits{i};
                % don't check existing fit names against itself
                if checkFitName(fit, targetFitdev) && strcmpi(fit.FitName, fitName)
                    taken = true;
                    break;
                end
            end
        end
        
        function sessionChanged(this)
            notify (this, 'FitsManagerSessionChanged'); 
        end
    end
     
    methods(Access = private )
        function this = addFit( this, hFitdev)
            % Add the given Fitdev to the list of fits
            this.Fits{end+1} = hFitdev;
            % Assign the name validator method
         
            hFitdev.CheckFitNameFcn = @(name, hFitdev)this.checkFitNameFcn(name, hFitdev);
            % Update the listeners.
            updateFitdevListeners( this )
        end

        function updateFitdevListeners( this )
            % Change the source of the Fitdev Listeners to match the list of
            % Fitdevs. If there are no listeners, then create new ones. If there
            % are no  fits, then delete the listeners.
            hFitdevs = this.Fits;
            
            if isempty( hFitdevs );
                % delete all listeners
                this.FitdevListeners = {};
                
            else
                % Need to create new listeners
                this.FitdevListeners = event.listener( hFitdevs, 'FitdevChanged', ...
                    @(s, e) this.sessionChanged() );
            end
        end
    end
end

function check = checkFitName(fit, targetFitdev)
    check = false;
    if isempty(targetFitdev)
        check = true;
    elseif (fit ~= targetFitdev)
        check = true;
    end
end

function [fit, i] = getFitFromName(Fits, fitName)
    fit = [];
    n = size(Fits, 2);
    for i = 1:n
        fit = Fits{i};
        if strcmp(fitName, fit.FitName)
            break;
        end
    end
end
 
function [fit, i] = getFitFromUUID(Fits, fitID)
    fit = [];
    n = size(Fits, 2);
    for i = 1:n
        fit = Fits{i};
        if fit.FitID.hashCode == fitID.hashCode
            break;
        end
    end
end
  
function [fitName, fitNameCount] = iGetUniqueFitName(this)
    fitName = sprintf('untitled fit %i', this.FitNameCount);
    fitNameCount = this.FitNameCount + 1;
    if isFitNameTaken(this, fitName, [])
        taken = true;
        while taken
            fitName = sprintf('untitled fit %i', fitNameCount);
            if isFitNameTaken(this, fitName, [])
                fitNameCount = fitNameCount + 1;
            else
                taken = false;
            end
        end
    end
end
  
function newName = iGetNameForDuplicatedFit(srcName, Fits)
    % Determine a new, unique name.
    name = srcName;
    taken = 1;
    i=1;

    ind=findstr(name,' copy ');
    if ~isempty(ind)
        name=name(1:(ind(end)-1));
    end

    % search for first unique name
    while taken
        newName = sprintf('%s copy %i',name,i);
        if iFitNameExists(newName, Fits)
            i = i+1;
        else
            taken = 0;
        end
    end
end
  
function nameExists = iFitNameExists(name, Fits)
    nameExists = false;
    n = size(Fits, 2);
    for i = 1:n
        if strcmp(Fits{i}.FitName, name)
            nameExists = true;
            break;
        end
    end
end
