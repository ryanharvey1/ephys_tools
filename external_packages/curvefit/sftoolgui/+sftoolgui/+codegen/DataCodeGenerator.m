classdef( Sealed ) DataCodeGenerator < handle
    % DataCodeGenerator   Generate code for data selected in SFTOOL
    
    %   Copyright 2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.1.2.1 $    $Date: 2011/07/18 00:31:04 $
    
    % The values of these properties are specified in the static methods
    % that wrap the constructor. They are used to specify the differences
    % in the generated code between fitting and validation data.
    properties( Access = private )
        XNameToken = 'undefined';
        YNameToken = 'undefined';
        ZNameToken = 'undefined';
        WNameToken = 'undefined';
        
        XVariableToken = 'undefined';
        YVariableToken = 'undefined';
        ZVariableToken = 'undefined';
        WVariableToken = 'undefined';
        
        XInputHelpFcn = [];
        YInputHelpFcn = [];
        YOutputHelpFcn = [];
        ZOutputHelpFcn = [];
        
        IgnoreWeights = false;
    end
    
    methods( Access = private )
        function cg = DataCodeGenerator()
            % DataCodeGenerator   Constructor
            %
            % This constructor is private. Please use one of the static
            % methods to create an instance of this class:
            %   sftoolgui.codegen.DataCodeGenerator.forFitting()
            %   sftoolgui.codegen.DataCodeGenerator.forValidation()
        end
    end
    
    % These static methods are the public interface to the constructor. In
    % addition to instantiating an object, they also need to specify the
    % values for the properties that distinguish fitting from validation.
    methods( Static )
        function cg = forFitting()
            % forFitting   Create an sftoolgui.codegen.DataCodeGenerator for
            % fitting data.
            cg = sftoolgui.codegen.DataCodeGenerator();
            
            cg.XNameToken = '<x-name>';
            cg.YNameToken = '<y-name>';
            cg.ZNameToken = '<z-name>';
            cg.WNameToken = '<w-name>';
            
            cg.XVariableToken = '<x-input>';
            cg.YVariableToken = '<y-input>';
            cg.ZVariableToken = '<z-output>';
            cg.WVariableToken = '<weights>';
            
            cg.XInputHelpFcn = @iXInputHelp;
            cg.YInputHelpFcn = @iYInputHelp;
            cg.YOutputHelpFcn = @iYOutputHelp;
            cg.ZOutputHelpFcn = @iZOutputHelp;
            
            cg.IgnoreWeights = false;
        end
        
        function cg = forValidation()
            % forValidation   Create an sftoolgui.codegen.DataCodeGenerator for
            % validation data.
            cg = sftoolgui.codegen.DataCodeGenerator();
            
            cg.XNameToken = '<xv-name>';
            cg.YNameToken = '<yv-name>';
            cg.ZNameToken = '<zv-name>';
            cg.WNameToken = 'DO NOT USE';
            
            cg.XVariableToken = '<validation-x>';
            cg.YVariableToken = '<validation-y>';
            cg.ZVariableToken = '<validation-z>';
            cg.WVariableToken = 'DO NOT USE';
            
            cg.XInputHelpFcn = @iXValidationHelp;
            cg.YInputHelpFcn = @iYValidationHelp;
            cg.YOutputHelpFcn = @iYValidationHelp;
            cg.ZOutputHelpFcn = @iZValidationHelp;
            
            cg.IgnoreWeights = true;
        end
    end
    
    methods
        function generateCode( cg, data, mcode )
            % generateCode   Generate code for an sftoolgui.Data
            %
            %   generateCode( CG, DATA, MCODE ) uses CG to generate code
            %   for DATA. The generated code is added to MCODE.
            %
            %   Types:
            %       CG     sftoolgui.codegen.DataCodeGenerator
            %       DATA   sftoolgui.Data
            %       MCODE  sftoolgui.codgen.MCode
            setFitVariables(   cg, data, mcode );
            addHelpComments(   cg, data, mcode );
            addPrepareCommand( cg, data, mcode );
        end
    end
    
    methods( Access = private )
        function tf = usingWeights( cg, data )
            % usingWeights   True for sftoolgui.Data that are using weights
            %
            %   This method takes into account the "IgnoreWeights" property
            %   as well as the weights in the sftoolgui.Data.
            [~, ~, ~, wName] = getNames( data );
            % Weights are being used if
            % ... the weights have a non-empty name
            % ... and the code generator has not been told to ignore them
            tf = ~isempty( wName ) && ~cg.IgnoreWeights;
        end
        
        function setFitVariables( cg, data, mcode )
            % setFitVariables   Set the names of the fit variables from the
            % sftoolgui.data into the sftoolgui.codgen.MCode
            [xName, yName, zName, wName] = getNames( data );
            % All valid data will have x and y variables
            iSetFitVariable( mcode, cg.XNameToken, xName );
            setFitVariable( mcode, cg.YNameToken, yName );
            % Curve data will not have a z variable
            iSetFitVariable( mcode, cg.ZNameToken, zName );
            % If we are using weights, then we need to register those
            if usingWeights( cg, data )
                setFitVariable( mcode, cg.WNameToken, wName );
            end
        end
        
        function addHelpComments( cg, data, mcode )
            % addHelpComments   Add the names of the variables from the
            % sftoolgui.data into help comments of the
            % sftoolgui.codgen.MCode
            [xName, yName, zName, wName] = getNames( data );
            
            % If there is an x variable, then it is always treated as an input
            if ~isempty( xName )
                addHelpComment( mcode, cg.XInputHelpFcn( xName ) );
            end
            
            % If we have curve data ...
            if isCurveDataSpecified( data )
                % .. then y should be the output
                addHelpComment( mcode, cg.YOutputHelpFcn( yName ) );
            else
                % ... but with surface data ...
                % ... y should be the input ...
                addHelpComment( mcode, cg.YInputHelpFcn( yName ) );
                % ... and z should be the output
                addHelpComment( mcode, cg.ZOutputHelpFcn( zName ) );
            end
            
            % If we are using weights add a comment for those
            if usingWeights( cg, data )
                addHelpComment( mcode, iWeightsHelp( wName ) );
            end
        end
        
        function addPrepareCommand( cg, data, mcode )
            % addHelpComments   Add a prepare command for the
            % sftoolgui.data into the sftoolgui.codgen.MCode
            
            xNameToken = xNameTokenForPrepareCommand( cg, data );
            
            if isCurveDataSpecified( data )
                if usingWeights( cg, data )
                    addFitCode( mcode, sprintf( ...
                        '[%s, %s, %s] = prepareCurveData( %s, %s, %s );', ...
                        cg.XVariableToken, cg.YVariableToken, cg.WVariableToken, ...
                        xNameToken,        cg.YNameToken,     cg.WNameToken ) );
                    
                else % ... not using weights
                    addFitCode( mcode, sprintf( ...
                        '[%s, %s] = prepareCurveData( %s, %s );', ...
                        cg.XVariableToken, cg.YVariableToken, ...
                        xNameToken,        cg.YNameToken ) );
                end
                
            else % ... surface data is specified
                if usingWeights( cg, data )
                    addFitCode( mcode, sprintf( ...
                        '[%s, %s, %s, %s] = prepareSurfaceData( %s, %s, %s, %s );', ...
                        cg.XVariableToken, cg.YVariableToken, cg.ZVariableToken, cg.WVariableToken, ...
                        xNameToken,        cg.YNameToken,     cg.ZNameToken,     cg.WNameToken ) );
                    
                else % ... not using weights
                    addFitCode( mcode, sprintf( ...
                        '[%s, %s, %s] = prepareSurfaceData( %s, %s, %s );', ...
                        cg.XVariableToken, cg.YVariableToken, cg.ZVariableToken, ...
                        xNameToken,        cg.YNameToken,     cg.ZNameToken ) );
                end
            end
        end
        
        function xNameToken = xNameTokenForPrepareCommand( cg, data )
            xValues = getValues( data );
            if isempty( xValues )
                xNameToken = '[]';
            else
                xNameToken = cg.XNameToken;
            end
        end
        
    end
end

% The following internal (sub) functions define the text that should go
% into the help comment of the generated code.
function str = iXInputHelp( name )
str = sprintf( '    X Input : %s', name );
end

function str = iYInputHelp( name )
str = sprintf( '    Y Input : %s', name );
end

function str = iYOutputHelp( name )
str = sprintf( '    Y Output: %s', name );
end

function str = iZOutputHelp( name )
str = sprintf( '    Z Output: %s', name );
end

function str = iWeightsHelp( name )
str = sprintf( '    Weights : %s', name );
end

function str = iXValidationHelp( name )
str = sprintf( '    Validation X: %s', name );
end

function str = iYValidationHelp( name )
str = sprintf( '    Validation Y: %s', name );
end

function str = iZValidationHelp( name )
str = sprintf( '    Validation Z: %s', name );
end

function iSetFitVariable( mcode, token, name )
% iSetFitVariable   Set a fit variable in code, but only if the NAME of the
% variable is not empty
if ~isempty( name )
    setFitVariable( mcode, token, name );
end

end

