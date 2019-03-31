classdef Fitdev < handle
    % Fitdev Surface Fitting Tool Fit Developer
    
    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.23.2.2 $    $Date: 2011/07/18 00:30:50 $ 

    events
        FittingDataUpdated
        ValidationDataUpdated
        FitFitted
        FitStarted
        FitNameUpdated
        CoefficientOptionsUpdated
        % FitTypeFitValuesUpdated is fired when a value that would require
        % a refit changes. These include:
        %   - data (both new data and new exclusions to current data)
        %   - fit options
        %   - fittype     
        FitTypeFitValuesUpdated
        % FitdevChanged is sent whenever Fitdev is changed to assist with
        % tracking session changes. It is expected that not much action
        % takes place in response to the FitdevChanged event as it is sent
        % frequently.
        FitdevChanged
        ExclusionsUpdated
        DimensionChanged
    end
    
    properties(SetAccess = 'private', GetAccess = 'private')
        % Version - class version number
        Version = 4;
        
        CurrentDefinition ;
        SurfaceDefinition ;
        CurveDefinition ;
        
        % ExclusionsStore - storage for exclusions. Empty implies all false.
        ExclusionsStore = [] ;
    end
    
    properties (Dependent)
        Exclusions ;
    end
    
    properties (Hidden, Dependent)
        % The following properties are here for compatibility and
        % convenience. They should not be used going forward. Instead, use
        % CurrentDefinition.Type and CurrentDefinition.Options or the
        % public method getCurrentDefinition which returns both properties.
        FitTypeObject ;
        FitOptions ;
    end
    
    properties(Transient)
        % CHECKFITNAMEFCN(NAME, THIS) should return empty if the NAME is a
        % valid fit name. Otherwise, it must return a struct with the
        % following fields: 'identifier', 'message'.
        CheckFitNameFcn ;
    end
    
    properties
        AutoFitEnabled = true;
        FitName ;
        FitID ;
        Fit ;
        Goodness ;
        ValidationGoodness ;
        Output ;
        FittingData ;
        ValidationData ;
        WarningStr ;
        ErrorStr ;
        ConvergenceStr ;
        FitState ;   % SFUtilities enum FitTypeObject: INCOMPLETE GOOD ERROR WARNING
        
        CustomEquation ;
    end
    
    properties (Dependent = true, SetAccess = private)
        ResultsTxt ;
        FitTypeString;
    end
    
    methods
        function this = Fitdev(name, id, sftoolguiData)
            % FITDEV Constructor for the Surface Fitting Tool Fit Developer
            %
            %   this = FITDEV(NAME, ID, SFTOOLGUIDATA)
            this.FitName = name;
            this.FitID = id;
            if isa(sftoolguiData, 'sftoolgui.Data')
                this.FittingData = sftoolguiData;
            elseif isempty(sftoolguiData)
                this.FittingData = sftoolgui.Data();
            else
                error(message('curvefit:Fitdev:InvalidInput'));
            end
            this.ValidationData = sftoolgui.Data();
            this.FitState = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.INCOMPLETE;
            
            % Create the default FitDefinition for surfaces.
            this.SurfaceDefinition = iCreateDefaultSurfaceDefinition();
            
            % Create the default FitDefinition for curves.
            this.CurveDefinition =  iCreateDefaultCurveDefinition();
            
            % Set current FitDefinition based on the data
            if isCurveDataSpecified( this.FittingData )
                this.CurrentDefinition = this.CurveDefinition;
            else
                % The default is the surface definition even if the data
                % specifies neither curve nor surface data.
                this.CurrentDefinition = this.SurfaceDefinition;
            end
            
            fit(this);
        end
        
        function set.FitName(this, name)
            if isempty(name)
                error(message('curvefit:sftoolgui:Fitdev:emptyFitName'));
            elseif ~ischar(name)
                error(message('curvefit:sftoolgui:Fitdev:NameNotString'));
            elseif ~isempty(this.CheckFitNameFcn) %#ok<MCSUP>
                error( this.CheckFitNameFcn (name, this)); %#ok<MCSUP>
            end
            % If we get here, we have a valid name
            this.FitName = name;
            notify (this, 'FitNameUpdated', sftoolgui.FitEventData( this));
        end
        
        function [type, options] = getCurveDefinition(this)
            % Get the CurveDefinition type and options.
            type = this.CurveDefinition.Type;
            options = this.CurveDefinition.Options;
        end
        
        function [type, options] = getSurfaceDefinition(this)
            % Get the SurfaceDefinition type and options
            type = this.SurfaceDefinition.Type;
            options = this.SurfaceDefinition.Options;
        end
        
        function [type, options] = getCurrentDefinition(this)
            % Get the CurrentDefinition type and options
            type = this.CurrentDefinition.Type;
            options = this.CurrentDefinition.Options;
        end
        
        function tf = isCurve(this)
            % Return true if Curve is the current definition, false
            % otherwise.
            tf = this.CurrentDefinition == this.CurveDefinition;
        end
        
        function ft = get.FitTypeString( this )
            if strcmpi(category(this.CurrentDefinition.Type), 'custom')
                ft = this.CustomEquation;
            else
                ft = type( this.CurrentDefinition.Type );
            end
        end
        
        function clearExclusions(this)        
            this.Exclusions = iInitializeExclusionArray(this.FittingData);
        end
        
        function toggleExclusion(this, index)
            this.Exclusions(index) = ~this.Exclusions(index);
        end
        
        function updateFittingData(this, evt)
            this.FittingData = setVariableFromBaseWS(this.FittingData, char(evt.getDataName), char(evt.getFieldName));
            doFittingDataUpdatedActions(this);
        end
        
        function clearFit(this)
            this.Fit = [];
            this.Goodness = [];
            this.ValidationGoodness = [];
            this.Output = [];
            this.WarningStr = [];
            this.ConvergenceStr = [];
            % If FitTypeObject is empty, ErrorStr can indicate a Fittype
            % error (generally due to an invalid custom equation). We want
            % to display this if users attempt a fit. Therefore, clear
            % ErrorStr only if FitTypeObject is not empty.
            if ~isempty(this.CurrentDefinition.Type)
                this.ErrorStr = [];
                this.FitState = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.INCOMPLETE;
            end
        end
        
        function updateValidationData(this, evt)
            this.ValidationData = setVariableFromBaseWS(this.ValidationData, char(evt.getDataName), char(evt.getFieldName));
            doValidationDataUpdatedActions(this);
        end
        
        function this = updateDataFromADifferentFit(this, srcFit)
            % Get the new data
            this.ValidationData = copy(srcFit.ValidationData);
            this.FittingData = copy(srcFit.FittingData);
            
            % doFittingDataUpdatedActions must be called before
            % doValidationDataUpdatedActions because
            % doValidationDataUpdatedActions needs information from the
            % fit.
            doFittingDataUpdatedActions(this);
            doValidationDataUpdatedActions(this);
        end
        
        function updateAutoFit(this, state)
            this.AutoFitEnabled = state;
            if this.AutoFitEnabled
                fit(this);
            end
            notify (this, 'FitdevChanged');
        end
        
        function updateCurveFittype(this, ft, fopts, customEquation, errorStr)
            % Update the curve fittype.
            this.CurveDefinition.Type = ft;
            this.CurveDefinition.Options = fopts;
            reconcileCoefficientOptions(this, this.CurveDefinition);
            this.CustomEquation = customEquation;
            this.ErrorStr = errorStr;
            
            if isCurve(this)
                updateFit(this);
            end
        end
        
        function updateCurveFitOptions(this, fopts)
            % Update the curve fit with new options.
            this.CurveDefinition.Options = fopts;
            if isCurve(this)
                updateFit(this);
            end
        end
        
        function updateSurfaceFittype(this, ft, fopts, customEquation, errorStr)
            % Update the surface fittype.
            this.SurfaceDefinition.Type = ft;
            this.SurfaceDefinition.Options = fopts;
            reconcileCoefficientOptions(this, this.SurfaceDefinition);
            this.CustomEquation = customEquation;
            this.ErrorStr = errorStr;
            
            if ~isCurve(this)
                updateFit(this);
            end
        end
        
        function updateSurfaceFitOptions(this, fopts)
            % Update the surface fit with new options.
            this.SurfaceDefinition.Options = fopts;
            if ~isCurve(this)
                updateFit(this);
            end
        end
        
        function saveFitToWorkspace(this)
            checkLabels = {'Save fit to MATLAB object named:', ...
                'Save goodness of fit to MATLAB struct named:', ...
                'Save fit output to MATLAB struct named:'};
            items = {this.Fit, this.Goodness, this.Output};
            varNames = {'fittedmodel', 'goodness', 'output'};
            export2wsdlg(checkLabels, varNames, items, 'Save Fit to MATLAB Workspace');
        end
        
        function tf = isFittingDataValid(this)
            % isFittingDataValid returns true if fitting data are specified
            % and data have compatible sizes. Otherwise it returns false.
            tf = iIsDataValid(this.FittingData);
        end
        
        function tf = isValidationDataValid(this)
            % isValidationDataValid   True for valid validation data
            %
            % isValidationDataValid returns true if validation data are
            % specified and have compatible sizes. Otherwise it returns
            % false.
            tf = iIsDataValid( this.ValidationData ) && isFittingValidationDataCompatible( this );
        end
        
        function this = fit(this)
            if isFittingDataValid(this)
                if ~isempty(this.CurrentDefinition.Type)
                    notify (this, 'FitStarted', sftoolgui.FitEventData( this));
                    
                    [input, output, w] = iGetInputOutputWeights(this.FittingData);
                    
                    % Assign weights and Exclusions into the fit options
                    fo = this.CurrentDefinition.Options;
                    
                    fo.Weights = w;
                    fo.Exclude = this.Exclusions;
                    
                    [this.Fit this.Goodness this.Output this.WarningStr this.ErrorStr this.ConvergenceStr] = ...
                        fit(input, output, this.CurrentDefinition.Type, fo );
                    
                    [this.ValidationGoodness.sse, this.ValidationGoodness.rmse] = getValidationGoodness( this );
                    if ~isempty(this.ErrorStr)
                        this.FitState = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.ERROR;
                    elseif ~isempty(this.WarningStr) || ~isempty(getMessageString(this.FittingData))
                        this.FitState = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.WARNING;
                    else
                        this.FitState = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.GOOD;
                    end
                else
                    % Should be an error, but if not, it is probably an invalid
                    % custom equation. Therefore set state to INCOMPLETE.
                    if ~isempty(this.ErrorStr)
                        this.FitState = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.ERROR;
                    else
                        this.FitState = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.INCOMPLETE;
                    end
                end
                
                notify (this, 'FitFitted', sftoolgui.FitEventData( this));
                notify (this, 'FitdevChanged');
            else
                this.Fit = [];
            end
        end
        
        function results = get.ResultsTxt( this )
            if isempty( this.Fit )
                % Empty results?
                results = '';
            else
                fitResults = genresults( ...
                    this.Fit, this.Goodness, this.Output, '', ...
                    this.ErrorStr, this.ConvergenceStr );
                
                validationString = makeValidationString( this );
                
                results = sprintf( '%s\n', fitResults{:}, validationString );
            end
        end
        
        function duplicatedFit = createADuplicate(this, newName)
            duplicatedFitUUID = javaMethodEDT('randomUUID', 'java.util.UUID');
            duplicatedFit = sftoolgui.Fitdev(newName, duplicatedFitUUID, sftoolgui.Data());
            duplicatedFit.AutoFitEnabled = this.AutoFitEnabled;
            duplicatedFit.Fit = this.Fit;
            duplicatedFit.Goodness = this.Goodness;
            duplicatedFit.ValidationGoodness = this.ValidationGoodness;
            duplicatedFit.Output = this.Output;
            duplicatedFit.FittingData = copy(this.FittingData);
            duplicatedFit.ValidationData = copy(this.ValidationData);
            duplicatedFit.WarningStr = this.WarningStr;
            duplicatedFit.ErrorStr = this.ErrorStr;
            duplicatedFit.ConvergenceStr = this.ConvergenceStr;
            duplicatedFit.FitState = this.FitState;
            duplicatedFit.ExclusionsStore = this.ExclusionsStore;
            duplicatedFit.CheckFitNameFcn = this.CheckFitNameFcn;
            duplicatedFit.CustomEquation = this.CustomEquation;
            copyDefinitions(this, duplicatedFit);
        end
        
        function generateMCode( this, mcode )
            % GENERATEMCODE   Generate Code for a Fitdev
            %
            %    GENERATEMCODE( H, CODE ) generates code for the given
            %    Fitdev, H, and adds it the code object CODE.
            INCOMPLETE = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.INCOMPLETE;
            
            if this.FitState == INCOMPLETE
                % Cannot generate code for "incomplete" fits
                addFitComment( mcode, sprintf( 'Cannot generate code for fit ''%s'' because the data selection is incomplete.', this.FitName ) );
                
            else % ... the fit is complete and we can generate code
                
                % Add code for the fitting data
                addHelpComment( mcode, sprintf( 'Data for ''%s'' fit:', this.FitName ) );
                iGenerateCodeForFittingData( this.FittingData, mcode );
                
                % Add code for fittype
                iGenerateCodeForFittype( this.CurrentDefinition.Type, mcode );
                
                % Add code for fit options
                [optionsFcn, extraFitArgsFcn] = iCodeForFitOptions( mcode, this.CurrentDefinition.Options );
                
                % Add code for weights
                if iHaveWeights( this.FittingData );
                    optionsFcn( 'Weights', '<weights>' );
                end
                
                % Add argument for center and scale, aka Normalize
                if strcmpi( this.CurrentDefinition.Options.Normalize, 'on' )
                    optionsFcn( 'Normalize', '''on''' );
                end
                
                % Add code to exclude points
                iGenerateCodeForExclusions( this.Exclusions, mcode, optionsFcn );
                
                % Add code to do fit
                % ... add a blank line to separate sections of code
                addBlankLine( mcode );
                % ... add a comment for the call to fit
                addFitComment( mcode, xlate( 'Fit model to data.' ) );
                % ... and the code to do the actual fit
                if isCurve( this )
                    addFitCode( mcode, '[<fo>, <gof>] = fit( <x-input>, <y-input>, <ft>, <opts> );' );
                else
                    addFitCode( mcode, sprintf( '[<fo>, <gof>] = fit( [<x-input>, <y-input>], <z-output>, <ft>%s );', extraFitArgsFcn() ) );
                end
                
                % Add code for validation data & validating the fit
                if isAnyDataSpecified( this.ValidationData )
                    iGenerateCodeToValidateFit( this, mcode );
                end
            end
        end
        
        function bIsFitted = isFitted(this)
            % isFitted returns true if a fitting operation was attempted and did not error.
            bIsFitted = this.FitState == com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.GOOD || ...
                this.FitState == com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.WARNING;  % good or warning fit
        end
        
        function fitTypeObject = get.FitTypeObject(this)
            % get fittype
            fitTypeObject = this.CurrentDefinition.Type;
        end
        
        function fitOptions = get.FitOptions(this)
            % get fit options
            fitOptions = this.CurrentDefinition.Options;
        end
        
        function set.Version(this, version)
            % set.Version was created so that load would create a struct
            % for objects whose version number is less than 4.
            if version >= 4
                this.Version = version;
            else
                error('curvefit:sftoolgui:Fitdev:IncompatibleVersion', ...
                    'Version number must be greater than 3');
            end
        end
        
        function exclusions = get.Exclusions(this)
            if isempty(this.ExclusionsStore)
                exclusions = iInitializeExclusionArray(this.FittingData);
            else
                exclusions = this.ExclusionsStore;
            end
        end
        
        function set.Exclusions(this, exclusions)
            this.ExclusionsStore = exclusions;
            
            notify (this, 'ExclusionsUpdated', sftoolgui.FitEventData( this));
            
            % Recalculate startpoints if it is an option and the fittype
            % has a startpt function.
            ft = this.CurrentDefinition.Type;
            % Check for empty startpt function here since we want to
            % preserve user's startpoints if there is no startpt function.
            % (calcStartPoints checks for presence of a startpt function
            % but will create random startpoints if the startpt function is
            % empty.)
            if isprop(this.CurrentDefinition.Options, 'StartPoint') && ~isempty(startpt(ft))
                this.CurrentDefinition.Options.StartPoint =  calcStartPoints(this, this.CurrentDefinition);
                notify (this, 'CoefficientOptionsUpdated', sftoolgui.FitEventData(this));
            end
            % Update fit
            updateFit(this);
        end
    end
    
    methods(Static)
        function obj = loadobj(this)
            
            % Get the version number.
            if isstruct(this) && ~isfield(this, 'Version')
                version = 0;
            else
                version = this.Version;
            end
            
            % Update THIS
            if version < 2
                this = iUpdateToV2(this);
            end
            
            if version < 3
                this = iUpdateV2ToV3(this);
            end
            
            if version < 4
                this = iUpdateV3ToV4(this);
            end
            
            % Create a new Fitdev object with THIS values
            try
                obj = iCreateLoadedObject(this);
            catch ignore %#ok<NASGU>
                obj = this;
                warning(message('curvefit:sftoolgui:Fitdev:unexpectedLoadedObject'));
            end
        end
    end
    
    methods(Access = private)
        function updateFit(this)
            % update Fit information that is relevant for both curves and
            % surface.
            % If there is an error ...
            if ~isempty(this.ErrorStr)
                % then update the FitState
                this.FitState = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.ERROR;
            end
            % call clearFit() after setting FitTypeObject because in
            % clearFit(), ErrorStr is conditionally cleared depending on
            % FitTypeObject state.
            clearFit(this);
            
            % Make sure FitTypeFitValuesUpdated event gets sent before
            % fitting
            notify (this, 'FitTypeFitValuesUpdated', ...
                sftoolgui.FitEventData(this));
            notify (this, 'FitdevChanged');
            if this.AutoFitEnabled
                fit(this);
            end
        end
        
        function doFittingDataUpdatedActions(this)
            % doFittingDataUpdatedActions does tasks required when fitting
            % data has changed.
            clearFit(this);
            previousDefinition = this.CurrentDefinition;
            % Update currentDefinition
            if isCurveDataSpecified(this.FittingData)
                this.CurrentDefinition = this.CurveDefinition;
            else
                this.CurrentDefinition = this.SurfaceDefinition;
            end
            % Dimension has changed
            if previousDefinition ~= this.CurrentDefinition
                notify (this, 'DimensionChanged');
            end          

            % Initialize the ExclusionsStore (Don't set this.Exclusions
            % here as that action does either some actions that we do not
            % want in this case or some actions that are redundant.)           
            this.ExclusionsStore = iInitializeExclusionArray(this.FittingData);
            
            reconcileCoefficientOptions(this, this.CurrentDefinition);
            
            notify (this, 'FittingDataUpdated', sftoolgui.FitEventData(this));
            notify (this, 'FitdevChanged');
            updateFit(this);
        end
        
        function doValidationDataUpdatedActions(this)
            % doValidationDataUpdateActions does tasks required when validation data has
            % changed.
            this.ValidationGoodness = [];
            if ~isempty(this.Fit)
                [this.ValidationGoodness.sse, this.ValidationGoodness.rmse] = getValidationGoodness( this );
            end
            notify (this, 'ValidationDataUpdated', sftoolgui.FitEventData(this));
            notify (this, 'FitdevChanged');
        end
               
        function startPoints =  calcStartPoints(this, fitDefinition)
            % calcStartPoints calculates startPoints for fittypes that have
            % startpt methods unless fitting data is invalid. Otherwise it
            % will return random startpoints.
            ft = fitDefinition.Type;
            spfun = startpt(ft);
            if isempty(spfun) || ~isFittingDataValid(this)
                startPoints = iCreateStartPoints(ft);
            else
                c = constants(ft);
                [inputs, output] = iGetInputOutputWeights(this.FittingData);
                % Remove excluded points from the data
                inputs = inputs(~this.Exclusions, :);
                output = output(~this.Exclusions);
                inputData = iCenterAndScaleInputs(inputs, fitDefinition.Options);
                % suppress the warning about power x needing to be positive
                warningState = warning('off', 'curvefit:fittype:sethandles:xMustBePositive');
                warningCleanup = onCleanup( @() warning( warningState) );
                startPoints = spfun(inputData, output, c{:});
            end
        end
        
        function reconcileCoefficientOptions(this, fitDefinition)
            % reconcileCoefficientOptions will ensure that startpoints,
            % upper bounds and lower bounds are all the same length. If the
            % fittype has a startpt method, it will use that to calculated
            % the startpoints. Otherwise it will just ensure that the
            % startpoints length is compatible with the number of
            % coefficients.
            
            ft = fitDefinition.Type;
            if isprop(fitDefinition.Options, 'StartPoint')
                fitDefinition.Options.StartPoint = calcStartPoints(this, fitDefinition);
            end
            if isprop(fitDefinition.Options, 'Lower')
                [lower, upper] = iCreateBounds(ft);
                fitDefinition.Options.Lower = lower;
                % If there are lower bounds, there will also be upper
                % bounds.
                fitDefinition.Options.Upper = upper;
            end
            % If there are lower bounds, there will also be upper bounds,
            % so check for just lower.
            if isprop(fitDefinition.Options, 'StartPoint') || isprop(fitDefinition.Options, 'Lower')
                notify (this, 'CoefficientOptionsUpdated', sftoolgui.FitEventData(this));
            end
        end
        
        function copyDefinitions(this, duplicatedFit)
            % copyDefinitions will create copies of surface and curve
            % FitDefinitions and assign the CurrentDefinition.
            % ('Type' is a fittype, which is a value property, so we do not
            % need to create a copy as we do with 'Options'.)
            duplicatedFit.SurfaceDefinition = sftoolgui.FitDefinition(this.SurfaceDefinition.Type, copy(this.SurfaceDefinition.Options));
            duplicatedFit.CurveDefinition = sftoolgui.FitDefinition(this.CurveDefinition.Type, copy(this.CurveDefinition.Options));
            if this.CurrentDefinition == this.CurveDefinition
                duplicatedFit.CurrentDefinition = duplicatedFit.CurveDefinition;
            else
                duplicatedFit.CurrentDefinition = duplicatedFit.SurfaceDefinition;
            end
        end
        
        function str = makeValidationString(this)
            % makeValidationString returns one of the following:
            % 1) validation statistics and possibly information about number of outside
            %    points (if validation goodness was calculated).
            % 2) an "incompatible data" message (if validation data is valid, but
            %    validation data and fitting data are incompatible).
            % 3) an empty string (if neither of the above cases applies).
            
            % Get the goodness of validation statistics.
            [sse, rmse, numOutside] = getValidationGoodness(this);
            
            % A non-empty sse indicates that validation goodness was calculated, which
            % means all of the following conditions were met:
            % 1) There was a fit object (which implies valid fitting data).
            % 2) Validation data was valid.
            % 3) Validation data and fitting data were compatible. (Data are
            %    incompatible if curve data is selected for fitting and surface data is
            %    selected for validation or surface data is selected for fitting and
            %    curve data is selected for validation).
            if ~isempty(sse)
                % Get message with statistics and possibly information about outside
                % points.
                if numOutside > 0
                    outsideMessage = sprintf( ...
                        '%d validation points outside domain of data\n', ...
                        numOutside );
                else
                    outsideMessage = '';
                end
                
                str = sprintf('SSE: %g\n  RMSE: %g\n%s', sse, rmse, outsideMessage );
                
                % An empty sse indicates validation goodness was not calculated. Are both
                % fitting data and validation data valid but incompatible? If so, get the
                % incompatible message.
            elseif  iIsDataValid(this.ValidationData) ...
                    && isFittingDataValid(this) ...
                    && ~isFittingValidationDataCompatible(this)
                
                str = iGetIncompatibleFittingValidationMessage();
                
                % Otherwise, validation goodness was not calculated for some other reason.
                % Return an empty string.
            else
                str = '';
            end
            % If there is a message, add a header.
            if ~isempty(str)
                str = sprintf('\nGoodness of validation:\n  %s', str);
            end
        end
        
        function [sse, rmse, numOutside] = getValidationGoodness( this )
            % getValidationGoodness returns validation entries and the number of
            % outside points.
            
            % Get the fit object.
            fo = this.Fit;
            
            % If the fit object is not empty and the validation data is valid ...
            if ~isempty( fo ) && isValidationDataValid(this) 
                % ... then calculate the validation goodness.
                [input, output] = iGetInputOutputWeights(this.ValidationData);
                [sse, rmse, numOutside] = iCalculateValidationGoodness( fo, input, output );
            else
                % ...otherwise, sse and rmse are empty, and the number of outside points is 0.
                sse = [];
                rmse = [];
                numOutside = 0;
            end
        end
        
        function tf = isFittingValidationDataCompatible(this)
            % isFittingValidationDataCompatible returns true if fitting and
            % validation data are either both curve or both surfaces.
            tf = isSurfaceData(this) || isCurveData(this);
        end
        
        function tf =  isSurfaceData(this)
            % isSurfaceData returns true if surface data is specified for both Fitting
            % and Validation and false otherwise.
            
            tf = isSurfaceDataSpecified(this.FittingData) && ...
                isSurfaceDataSpecified(this.ValidationData);
        end
        
        function tf =  isCurveData(this)
            % isCurveData returns true if curve data is specified for both Fitting
            % and Validation and false otherwise.
            tf = isCurveDataSpecified(this.FittingData) && ...
                isCurveDataSpecified(this.ValidationData);
        end
    end
end

function this = iUpdateToV2(this)
% iUpdateToV2 adds or updates the 'CustomEquation' property. The
% 'CustomEquation' field was added for the 9a release (i.e. was not in the
% Prerelease version) but the 'Version' number was not increased. THIS is
% a struct that represents Fitdevs whose 'Version' property is < 2.

% If there isn't a 'CustomEquation' field
if ~isfield(this, 'CustomEquation')
    % ... then we need to add one
    this.CustomEquation = '';
end

% If this Fitdev is for a custom equation, but the 'CustomEquation' field
% is empty
if strcmp( category( this.FitTypeObject ), 'custom' ) && isempty( this.CustomEquation )
    % ... then we need to guess that the formula is close enough to what
    % the user typed for the custom equation (it will be accurate to
    % within some amount of whitespace)
    this.CustomEquation = formula( this.FitTypeObject );
end
end

function this = iUpdateV2ToV3(this)
% iUpdateV2ToV3 handles incompatible exclusions. If the exclusion array is
% longer than the longest data array, it probably means that the Data obj
% is V1 (before we removed NaNs and Infs) so in that case reinitialize
% array and warn the user if there were any exclusions.
[x, y, z] = getValues(this.FittingData);
maxDataLength = max([length(x), length(y), length(z)]);
excludeLength = length(this.Exclusions);
if excludeLength > maxDataLength
    if any(this.Exclusions)
        warning(message('curvefit:sftoolgui:Fitdev:ExclusionsLost', this.FitName));
    end
    this.Exclusions = [];
end
end

function this = iUpdateV3ToV4(this)
% iUpdateV3ToV4 creates a SurfaceDefinition from Fitdev Version 3
% FitTypeObject and FitOptions, creates the default CurveDefinition and
% sets the Current Definition to be Surface. It also sets the
% ExclusionsStore property from the old Exclusions property.
if isfield(this, 'FitTypeObject') && isfield(this, 'FitOptions')
    this.SurfaceDefinition = sftoolgui.FitDefinition(this.FitTypeObject, this.FitOptions);
else
    this.SurfaceDefinition = iCreateDefaultSurfaceDefinition();
end
this.CurveDefinition =  iCreateDefaultCurveDefinition();
this.CurrentDefinition = this.SurfaceDefinition;

if isfield(this, 'Exclusions')
    this.ExclusionsStore = this.Exclusions;
end
end

function obj = iCreateLoadedObject(this)
% iCreateLoadedObject creates an sftoolgui.Fitdev and assigns values of
% THIS to its properties. THIS is either a struct or an sftoolgui.Fitdev.
obj = sftoolgui.Fitdev(this.FitName, this.FitID, []);
obj.SurfaceDefinition = this.SurfaceDefinition;
obj.CurveDefinition = this.CurveDefinition;
obj.CurrentDefinition = this.CurrentDefinition;
obj.AutoFitEnabled = this.AutoFitEnabled;
obj.Fit = this.Fit;
obj.Goodness = this.Goodness;
obj.ValidationGoodness = this.ValidationGoodness;
obj.Output = this.Output;
obj.FittingData = this.FittingData;
obj.ValidationData = this.ValidationData;
obj.WarningStr = this.WarningStr;
obj.ErrorStr = this.ErrorStr;
obj.ConvergenceStr = this.ConvergenceStr;
obj.FitState = this.FitState;
obj.ExclusionsStore = this.ExclusionsStore;
obj.CustomEquation = this.CustomEquation;
end

function exclusions = iInitializeExclusionArray(fittingData)
%iInitializeExclusionArray creates an exclusions array. The length of the
%array is the maximum length of x, y and z.
if isCurveDataSpecified(fittingData)
    [x, y] = getCurveValues(fittingData);
    z = [];
else
    [x, y, z] = getValues(fittingData);
end
excludeLength = max([length(x), length(y), length(z)]);
exclusions = false(1, excludeLength);
end

function isValid = iIsDataValid(data)
isValid = (isSurfaceDataSpecified(data) || isCurveDataSpecified(data)) ...
    && areNumSpecifiedElementsEqual(data);
end

function tf = iHaveWeights( hData )
% iHaveWeights -- Determine if a Data object has weights
[~, ~, ~, wName] = getNames( hData );
tf = ~isempty(wName);
end

%--- Code Generation Sub-Functions
function [optionsFcn, extraFitArgsFcn] = iCodeForFitOptions( mcode, fopts )
% Generate code for fit options
%
% This function returns two function handles:
%
%     optionsFcn = @(p, v)   is a function for adding code for a fit option
%         parameter and value.
%
%     extraFitArgsFcn = @()  is a function that returns a string with a list of
%         addition arguments to pass the fit function.
%
% The optionsFcn will vary depending on whether a fit options object is required
% in the generated code.
requiresOptions = false;
optionsCode = '<opts> = fitoptions( <ft> );';

% Work out the default options
fn = fieldnames( fopts );
defaultOptions = fitoptions( 'Method', fopts.Method );
fn = setdiff( fn, 'Method' );

% Work out if any interesting options differ from the defaults and if so, what
% code should be used to generated them
for i = 1:length( fn )
    if isequal( defaultOptions.(fn{i}), fopts.(fn{i}) );
        % no need for code
    elseif ismember( fn{i}, {'Exclude', 'Normalize'} )
        % Code for Exclude and Normalize is handled elsewhere
    elseif isequal( fn{i}, 'Weights' )
        % The vector of weights needs to be assigned from user data and so
        % is handled elsewhere, however it also needs an explicit fit
        % options object to be used.
        requiresOptions = true;
    else
        requiresOptions = true;
        optionsCode = sprintf( '%s\n<opts>.%s = %s;', optionsCode, fn{i}, ...
            mat2str( fopts.(fn{i}) ) );
    end
end

if requiresOptions
    % If an options object is required, then return a handle that will add code
    % for that object
    optionsFcn = @nAddOptions;
    % Also need to declare a variable in the code to hold those options
    addVariable( mcode, '<opts>', 'opts' );
    % Add the code to set non-default values of the options
    addFitCode( mcode, optionsCode );
    % Set the extra fit arguments to include the fit options we have created
    extraFitArguments = ', <opts>';
else
    % If a fit options object is no required, then return a handle to a function
    % that add code to the fit arguments
    optionsFcn = @nAddFitArguments;
    % Initialize those arguments to empty.
    extraFitArguments = '';
end

extraFitArgsFcn = @nReturnExtraFitArguments;

% Add code to add an option to the fit-options object
    function nAddOptions( parameter, value )
        addFitCode( mcode, sprintf( '<opts>.%s = %s;', parameter, value ) );
    end

% Generate p-v pairs for the fit command
    function nAddFitArguments( parameter, value )
        extraFitArguments = sprintf( '%s, ''%s'', %s', extraFitArguments, parameter, value );
    end

% Return a string of arguments to give to the fit command
    function args = nReturnExtraFitArguments()
        args = extraFitArguments;
    end

end

function iGenerateCodeForFittingData( hData, mcode )
% iGenerateCodeForFittingData -- Generate Code for Fitting Data
theGenerator = sftoolgui.codegen.DataCodeGenerator.forFitting();
generateCode( theGenerator, hData, mcode );
end

function iGenerateCodeForFittype( ft, mcode )
% iGenerateCodeForFittype -- Generate code for a fittype
addVariable( mcode, '<ft>', 'ft' );

% Add a blank line to separate sections of code
addBlankLine( mcode );
% Add a comment for the fittype section
addFitComment( mcode, xlate( 'Set up fittype and options.' ) );

% Add the code for the fittype
if strcmpi( category( ft ), 'custom' )
    iAddCustomFittypeCode( ft, mcode );
    
elseif ismember( type( ft ), {'linearinterp', 'nearestinterp', 'cubicinterp'} )
    addFitCode( mcode, sprintf( '<ft> = ''%s'';', type( ft ) ) );
    
else
    addFitCode( mcode, sprintf( '<ft> = fittype( ''%s'' );', type( ft ) ) );
end
end

function iAddCustomFittypeCode( ft, mcode )
% iAddCustomFittypeCode   Add code for a custom equation fittype to the given
% MCode object.
independentVariables = indepnames( ft );
dependentVariables = dependnames( ft );

% The line for fittype is made up of four pieces:
% ... 1. the start (including the equation)
start = sprintf( '<ft> = fittype( ''%s'', ', formula( ft ) );

% ...2. independent variables
if length( independentVariables ) == 1
    independentCode = sprintf( '''independent'', ''%s'', ', independentVariables{1} );
    
else % length( indepVars ) == 2
    independentCode = sprintf( '''independent'', {''%s'', ''%s''}, ', independentVariables{:} );
end

% ... 3. Dependent variable (always exactly one)
dependentCode = sprintf( '''dependent'', ''%s''', dependentVariables{1} );

% .. 4. The end of the line is closing brace and semi-colon
finish = ' );';

% Concatenate the four pieces and add to the MCode
addFitCode( mcode, [start, independentCode, dependentCode, finish] );

end

function iGenerateCodeForExclusions( exclusions, mcode, optionsFcn )
% iGenerateCodeForExclusions -- Generate code to define an exclusion vector
% and assign it as one of the options
%
% Inputs:
%   exclusions: logical array where true indicates a point to exclude
%   mcode: the sftoolgui.codegen.MCode object to add code to
%   optionsFcn: a handle to a function to add options code to the MCode
%       object. This is one of the return arguments from
%       iCodeForFitOptions.
if any( exclusions )
    % Need to store the exclusion variable for use with plots
    addVariable( mcode, '<ex>', 'ex' );
    ex = find( exclusions );
    addFitCode( mcode, sprintf( '<ex> = excludedata( <x-input>, <y-input>, ''Indices'', %s );', mat2str( ex ) ) );
    optionsFcn( 'Exclude', '<ex>' );
end
end

function iGenerateCodeToValidateFit( this, mcode )
% iGenerateCodeToValidateFit   Generate the code required to validate a fit

% Add a blank line to separate sections of code
addBlankLine( mcode );

if isValidationDataValid( this )
    % If the validation data is valid ...
    % ... add a comment for the validation section
    addFitComment( mcode, xlate( 'Compare against validation data.' ) );
    % ... add code for the validation data
    iGenerateCodeForValidationData( this.ValidationData, mcode );
    % ... add code for the goodness of validation computation and display
    iGenerateCodeForGoodnessOfValidation( this.FitName, mcode, isCurveDataSpecified( this.ValidationData ) );

else
    % If the validation data is invalid ...                       
    % ... insert a comment telling the user that we couldn't generate code for them.
    addFitComment( mcode, sprintf( 'Cannot generate code for validating ''%s'' because the validation data is incomplete or incompatible with the fitting data.', this.FitName ) );
end
end

function iGenerateCodeForValidationData( hData, mcode )
% iGenerateCodeForValidationData -- Generate Code for Validation Data
theGenerator = sftoolgui.codegen.DataCodeGenerator.forValidation();
generateCode( theGenerator, hData, mcode );
end

function iGenerateCodeForGoodnessOfValidation( fitName, mcode, isCurve )
% iGenerateCodeForGoodnessOfValidation -- Generate the code that computes
% the goodness of validation and then displays it.

addVariable( mcode, '<residual>', 'residual' );
addVariable( mcode, '<sse>', 'sse' );
addVariable( mcode, '<rmse>', 'rmse' );
addVariable( mcode, '<nNaN>', 'nNaN' );

if isCurve
    addFitCode( mcode, '<residual> = <validation-y> - <fo>( <validation-x> );' );
else
    addFitCode( mcode, '<residual> = <validation-z> - <fo>( <validation-x>, <validation-y> );' );
end

addFitCode( mcode, '<nNaN> = nnz( isnan( <residual> ) );' );
addFitCode( mcode, '<residual>(isnan( <residual> )) = [];' );
addFitCode( mcode, '<sse> = norm( <residual> )^2;' );
addFitCode( mcode, '<rmse> = sqrt( <sse>/length( <residual> ) );' )

% Need to separate string construction from code
% construction for proper translation.
gofString  = xlate( 'Goodness-of-validation for ''''%s'''' fit' );
sseString  = xlate( 'SSE' );
rmseString = xlate( 'RMSE' );
outsideString = xlate( '%d points outside domain of data.' );

addFitCode( mcode, sprintf( 'fprintf( ''%s:\\n'', ''%s'' );', gofString, fitName ) );
addFitCode( mcode, sprintf( 'fprintf( ''    %s : %%f\\n'', <sse> );',  sseString  ) );
addFitCode( mcode, sprintf( 'fprintf( ''    %s : %%f\\n'', <rmse> );', rmseString ) );
addFitCode( mcode, sprintf( 'fprintf( ''    %s\\n'', <nNaN> );', outsideString ) );
end
%--- End of Code Generation Sub-Functions

function [sse, rmse, numOutside] = iCalculateValidationGoodness( fo, input, output )
% iCalculateValidationGoodness returns validation entries and number of
% outside points. It assumes that fo, input and output are all valid.
outputHat = fo(input);
% For some model types, esp. GRIDDATA, outputHat will be NaN for some
% points because they are outside the region of definition of the
% surface. We remove those points from the stats and warn the user.
outputHatIsNan = isnan( outputHat );
numOutside = nnz( outputHatIsNan );
if numOutside > 0
    output = output(~outputHatIsNan);
    outputHat = outputHat(~outputHatIsNan);
end
residual = output(:) - outputHat;
sse = norm( residual )^2;
rmse = sqrt( sse/length( outputHat ) );
end

function message = iGetIncompatibleFittingValidationMessage()
% iGetIncompatibleFittingValidationMessage returns an "incompatible fitting
% and validation data" message.
message = xlate('Not calculated. Incompatible fitting and validation data.');
end

function fitDefinition = iCreateDefaultCurveDefinition()
% iCreateDefaultCurveDefinition creates the default sftoolgui.FitDefinition
% for curves
defaultFittype = fittype('poly1');
defaultFitOptions = fitoptions(defaultFittype);
fitDefinition = sftoolgui.FitDefinition(defaultFittype, defaultFitOptions);
end

function fitDefinition = iCreateDefaultSurfaceDefinition()
% iCreateDefaultSurfaceDefinition creates the default
% sftoolgui.FitDefinition for surfaces and sets the 'Normalize'
% option to be 'on'.
defaultFittype = fittype( 'linearinterp', 'numindep', 2 );
defaultFitOptions = fitoptions(defaultFittype);
defaultFitOptions.Normalize = 'on';
fitDefinition = sftoolgui.FitDefinition(defaultFittype, defaultFitOptions);
end

function [lower, upper] = iCreateBounds(fitType)
% iCreateBounds creates arrays of lower and upper bounds. The length of the
% arrays is the same as the number of coefficients.
numcoeff = numcoeffs(fitType);
% If the fittype has specific lower bounds, use those.
fo = fitoptions(fitType);
if isprop(fo, 'Lower') && ~isempty(fo.Lower)
    lower = fo.Lower;
else
    lower = -Inf(1, numcoeff);
end
if isprop(fo, 'Upper') && ~isempty(fo.Upper)
    upper = fo.Upper;
else
    upper = Inf(1, numcoeff);
end
end

function startPoints = iCreateStartPoints(fitType)
% iCreateStartPoints creates an array of (random) startpoints that is the
% same length as the number of coefficients.
numcoeff = numcoeffs(fitType);
startPoints = rand(1, numcoeff);
end

function inputs = iCenterAndScaleInputs(inputs, fopts)
% iCenterAndScaleInputs will center and scale inputs if "Normalize" option
% is on.
if isprop(fopts, 'Normalize') && strcmpi(fopts.Normalize, 'on')
    inputs = curvefit.normalize(inputs);
end
end

function [input, output, weights] = iGetInputOutputWeights(data)
% getInputOutputWeights returns input, output and weights
% for both curve and non-curve data, which is assumed to have been tested
% for validity before this function is called.

if isCurveDataSpecified(data)
    [input, output, weights] = getCurveValues(data);
else
    [x, y, output, weights] = getValues(data);
    % concatenate x and y
    input = [x, y];
end
end
