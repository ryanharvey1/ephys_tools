classdef MCode < handle
    %MCODE   The MCode object represents code generated from a SFTOOL session
    %
    %   OBJ = MCODE
    
    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.8 $    $Date: 2011/05/09 00:39:58 $

    properties(GetAccess = 'public', Constant = true)
        % FitVariableTokens -- Tokens that can be used in code but that will be
        % replaced by the corresponding variable for the fit.
        FitVariableTokens = {'<x-name>', '<y-name>', '<z-name>', '<w-name>', '<xv-name>', '<yv-name>', '<zv-name>'};
    end
    properties(SetAccess = 'private', GetAccess = 'private')
        % FitResult -- variable name for the fit object(s) that get returned by the
        % generated code
        FitResult = 'fitresult';
        
        % GOFVariable -- variable name for the goodness-of-fit (GOF) structure
        % that gets returned by the generated code
        GOFVariable = 'gof';
        
        % Inputs -- List of inputs to the generated code
        Inputs = {};
        
        % FitVariables -- Names of the variables for each fit
        FitVariables;
        
        % HelpComments -- Text of the comment "help" block of the code
        HelpComments = '';
        
        % FitBlocks -- Cell-array of cell-strings.
        %
        %   FitBlocks{n} is a cell-array of strings and represents the code that
        %   fits and plots the n-th fit in the session.
        FitBlocks = {}
        
        % VariableTokens -- Cell-array of strings
        VariableTokens = {}
        % Variables -- Cell-array of strings
        Variables = {}
    end
    
    methods
        function h = MCode
            h.FitVariables = cell( 0, length( h.FitVariableTokens ) );
            
            % Declare all the variables that will be used for fitting
            addVariable( h, '<x-input>', 'xData' );
            addVariable( h, '<y-input>', 'yData' );
            addVariable( h, '<z-output>', 'zData' );
            addVariable( h, '<weights>', 'weights' );
            % ... and validation
            addVariable( h, '<validation-x>', 'xValidation' );
            addVariable( h, '<validation-y>', 'yValidation' );
            addVariable( h, '<validation-z>', 'zValidation' );
        end
        
        function startNewFit( h )
            % STARTNEWFIT   Notify MCode object of start of new fit
            h.FitBlocks{end+1} = {};
            h.FitVariables(end+1,:) = h.FitVariableTokens;
        end
        
        function setFitVariable( h, variable, name )
            % SETFITVARIABLE   Set a variable for this fit
            %
            %   setFitVariable( h, '<x-name>', name )
            %   setFitVariable( h, '<y-name>', name )
            %   setFitVariable( h, '<z-name>', name )
            %   setFitVariable( h, '<weights>', name )
            [tf, loc] = ismember( variable,  h.FitVariableTokens );
            if tf
                h.FitVariables{end,loc} = name;
                addInput( h, name );
            else
                error(message('curvefit:sftoolgui:MCode:InvalidFitVariable'));
            end
        end
        
        function addVariable( h, token, name )
            % ADDVARIABLE   Register a variable to the generated code
            if ismember( token, h.VariableTokens )
                % do nothing
            else
                h.VariableTokens{end+1} = token;
                h.Variables{end+1}      = name;
            end
        end
        
        function addInput( h, input )
            % ADDINPUT   Register an input to the generated code
            if ismember( input, h.Inputs )
                % do nothing
            else
                h.Inputs{end+1} = input;
            end
        end
        
        function addHelpComment( h, comment )
            % ADDHELPCOMMENT   Register a line in the help comments
            %
            %   See also: addFitComment
            h.HelpComments = sprintf( '%s\n%s', h.HelpComments, comment );
        end
        
        function addFitComment( h, comment )
            % ADDFITCOMMENT  Add a comment line to the code for fitting a surface
            %
            %   ADDFITCOMMENT( H, COMMENT ) adds the comment, COMMENT, to the
            %   MATLAB Code to be generated, H. COMMENT should be a char array
            %   representing a single line, i.e., it should not a '\n' in it.
            %
            %   A comment character followed by a space, '% ', is prepended to
            %   the comment before it is written to the file.
            %
            %   The COMMENT should be translated before passing into this
            %   method, e.g.,
            %       addFitComment( mcode, xlate( 'This is my comment' ) )
            %
            %   See also: addHelpComment
            h.FitBlocks{end}{end+1} = sprintf( '%% %s', comment );
        end
        
        function addCellHeader( h, line )
            % ADDCELLHEADER  Add a cell header to the code for fitting a surface 
            %
            %   ADDCELLHEADER( H, COMMENT ) adds the line, LINE, to the MATLAB
            %   Code to be generated, H, and makes it into a cell header. LINE
            %   should be a char array representing a single line, i.e., it
            %   should not a '\n' in it.
            %
            %   The cell-header markup followed by a space, '%% ', is prepended
            %   to the line before it is written to the file.
            h.FitBlocks{end}{end+1} = iMakeCellHeader( line );
        end
        
        function addBlankLine( h )
            % ADDBLANKLINE  Add a blank line into code
            %
            %  ADDBLANKLINE( H ) adds an empty line in the MATLAB code to be 
            %  generated, H.
            h.FitBlocks{end}{end+1} = '';
        end
        
        function addFitCode( h, code )
            % ADDFITCODE   Add a code for fitting a surface to MCode object
            %
            %   ADDFITCODE( H, CODE ) adds the code, CODE, to the MATLAB Code to
            %   be generated H. CODE should be a char array representing one of
            %   more lines. Separate multiple lines with a new line character,
            %   i.e., sprintf( '\n' ).
            %
            %   The code can include various tokens. In the generated code,
            %   these tokens will be replaced by appropriate variables name. The
            %   allowed tokens and their meanings are:
            %
            %       <fo>    fit object
            %       <gof>   goodness-of-fit
            %       <x-name> name of the x variable
            %       <y-name> name of the y variable
            %       <z-name> name of the z variable
            %       <w>     weights
            %
            %   Note that in the generated code, each of these variables may
            %   vary with each fit.
            h.FitBlocks{end}{end+1} = code;
        end
        
        function hFunction = coderoutine( h )
            % coderoutine   Convert an MCode object to a codegen.coderoutine object
            nFits = length( h.FitBlocks );
            
            % Set the names of the output and temporary variables
            finalizeVariables( h );
            
            % Create the coderoutine object
            hFunction = codegen.coderoutine;

            % Add header information
            % -- name the function.
            if nFits == 1
                hFunction.Name = 'createFit';
            else
                hFunction.Name = 'createFits';
            end
            % -- register input arguments
            for i = 1:length( h.Inputs )
                hFunction.addArgin( h.Inputs{i} );
            end
            % -- register output arguments
            hFunction.addArgout( h.FitResult );
            hFunction.addArgout( h.GOFVariable );
            % -- add help comments at top of file
            hFunction.Comment = generateHelpComment( h );
            hFunction.SeeAlsoList = {'fit', 'cfit', 'sfit'};

            % Add code to initialize output arguments
            if nFits > 1
            hFunction.addText( iMakeCellHeader( xlate( 'Initialization.' ) ) );
            hFunction.addText( '' ); % blank line
                hFunction.addText( iMakeComment( xlate( 'Initialize arrays to store fits and goodness-of-fit.' ) ) );
                hFunction.addText( sprintf( '%s = cell( %d, 1 );', h.FitResult, nFits ) );
                hFunction.addText( sprintf( '%s = struct( ''sse'', cell( %d, 1 ), ...', h.GOFVariable, nFits ) );
                hFunction.addText( '''rsquare'', [], ''dfe'', [], ''adjrsquare'', [], ''rmse'', [] );' );
            end
            
            % Add code for fitting and plotting
            % -- get names of output variables as they change for each fit
            if nFits == 1
                fitObjectName = @(i) h.FitResult;
                gofName       = @(i) h.GOFVariable;
            else
                fitObjectName = @(i) sprintf( '%s{%d}', h.FitResult, i );
                gofName       = @(i) sprintf( '%s(%d)', h.GOFVariable, i );
            end
            % -- add fitting and plotting code for each fit in turn
            for i = 1:nFits
                hFunction.addText( '' );
                for j = 1:length( h.FitBlocks{i} )
                    thisLine = h.FitBlocks{i}{j};
                    thisLine = strrep( thisLine, '<fo>', fitObjectName( i ) );
                    thisLine = strrep( thisLine, '<gof>', gofName( i ) );
                    thisLine = replaceTokens( h, thisLine, i );
                    hFunction.addText( thisLine );
                end
            end
            hFunction.addText( '' );
        end
    end
    
    methods(Access = 'private' )
        function finalizeVariables( h )
            % This method ensures that the variable names do not clash with any
            % of the inputs or any other variables that have been registered
            allVariables = h.Inputs;
            
            % FitResult
            h.FitResult = genvarname( h.FitResult, allVariables );
            allVariables{end+1} = h.FitResult;

            % Goodness of fit
            h.GOFVariable = genvarname( h.GOFVariable, allVariables );
            allVariables{end+1} = h.GOFVariable;
            
            % Temporary Variables
            h.Variables = genvarname( h.Variables, allVariables );
        end
        
        function line = replaceTokens( h, line, fitIndex )
            % Replace tokens in a line of code
            for i = 1:length( h.FitVariableTokens )
                line = strrep( line, h.FitVariableTokens{i}, h.FitVariables{fitIndex,i} );
            end
            for i = 1:length( h.Variables )
                line = strrep( line, h.VariableTokens{i}, h.Variables{i} );
            end
        end
        
        function comments = generateHelpComment( h )
            % generateHelpComment -- Generate the help comment that gets
            % inserted at the top of file.
            nFits = length( h.FitBlocks );
            if nFits == 1
                line1 = xlate( 'Create a fit.' );
                output1 = xlate( 'Output:' ); 
                output2 = sprintf( '    %s : a fit object representing the fit.', h.FitResult );
                output3 = sprintf( '    %s : structure with goodness-of fit info.', h.GOFVariable );
            else
                line1 = xlate( 'Create fits.' );
                output1 = xlate( 'Output:' );
                output2 = sprintf( '    %s : a cell-array of fit objects representing the fits.', h.FitResult );
                output3 = sprintf( '    %s : structure array with goodness-of fit info.', h.GOFVariable );
            end
            comments = sprintf( '%s\n%s\n%s\n%s\n%s', line1, h.HelpComments, output1, output2, output3 );
        end
    end
end

function comment = iMakeComment( text )
% iMakeComment -- Make a comment from a piece of text
%
% The purpose of this routine is to make translation easier. It may be called as
% follows: 
%   addText( hFunction, iMakeComment( xlate( 'my comment' ) ) )
% where hFunction us a codegen.coderoutine.
%
% The comment character plus a space, '% ' will be added to the TEXT to make the
% comment
comment = sprintf( '%% %s', text );
end

function cellHeader = iMakeCellHeader( text )
% iMakeComment -- Make a cell header from a piece of text
%
% The cell break symbol plus a space, '%% ' will be added to the TEXT to make
% the cell header.
cellHeader = sprintf( '%%%% %s', text );
end
