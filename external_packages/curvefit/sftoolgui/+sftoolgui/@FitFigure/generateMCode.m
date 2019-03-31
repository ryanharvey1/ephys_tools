function generateMCode( this, mcode )
%generateMCode Generate code for a Fit Figure
%
%   generateMCode( H, CODE ) generates code for the given fit
%   figure, H, and adds it the code object CODE.

%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.6 $    $Date: 2011/05/09 00:40:14 $

INCOMPLETE_FIT = com.mathworks.toolbox.curvefit.surfacefitting.SFFitState.INCOMPLETE;

startNewFit( mcode );
addCellHeader( mcode, sprintf( 'Fit: ''%s''.', this.HFitdev.FitName ) )
generateMCode( this.HFitdev, mcode );

if this.HFitdev.FitState ~= INCOMPLETE_FIT
    
    % Code to create a new figure
    codeForFigure = sprintf( 'figure( ''Name'', ''%s'' );', ...
        sftoolgui.codegen.stringLiteral( this.HFitdev.FitName ) );
    
    % Which panels do we need to generate code for?
    [nPlots, panel, comment] = iWhichPanels( ...
        this.HSurfacePanel, this.HResidualsPanel, this.HContourPanel);
    
    % Code to handle the layout of the given number of plots that will have code generated.
    switch nPlots
        case 0
            layout = cell(0);
        case 1
            layout{1} = codeForFigure;
        case 2
            addBlankLine( mcode );
            addFitComment( mcode, xlate( 'Create a figure for the plots.' ) );
            addFitCode( mcode, codeForFigure );
            layout{1} = 'subplot( 2, 1, 1 );' ;
            layout{2} = 'subplot( 2, 1, 2 );' ;
        case 3
            addBlankLine( mcode );
            addFitComment( mcode, xlate( 'Create a figure for the plots.' ) );
            addFitCode( mcode, codeForFigure );
            layout{1} = 'subplot( 2, 2, 2 );' ;
            layout{2} = 'subplot( 2, 2, 4 );' ;
            layout{3} = 'subplot( 1, 2, 1 );' ;
        otherwise
            % Should not be able to get to here.
            warning(message('curvefit:sftoolgui:FitFigure:InvalidState'));
            layout = cell( 1, nPlots );
    end
    
    % If there is validation data, then the generated code must work
    % out the axes limits
    iGenerateAxesLimitCode( this.HFitdev, mcode );
    
    % Generate code or each of the plots that can have code generated
    for i = 1:nPlots
        addBlankLine( mcode );
        addFitComment( mcode, comment{i} );
        addFitCode( mcode, layout{i} );
        generateMCode( panel{i}, mcode );
    end
    
    % Add a comment for the contour panel
    iGenerateContourForCurvesComment(this.HContourPanel, mcode);
end
end

function [numPlots, codeGenPanels, comments] = iWhichPanels( ...
    surfacePanel, residualsPanel, contourPanel)
% iWhichPanels   Work out which panels to generated code for.
%
% This function works out how many plots (numPlots) need to have code
% generated and which plots those are (codeGenPanels). It also returns the
% appropriate comment to introduce the generated code for each panel.

% For each panel where we want to generate code, we will increment this
% counter by one.
numPlots = 0;
% Just in case there are no panels to generate code for we will declare
% empty codeGenPanels and comments cell arrays.
codeGenPanels = {};
comments = {};

% For each panel, if we can generate code, increment the counter and append
% the panel to codeGenPanels and add an appropriate comment to comments.
if canGenerateCodeForPlot(surfacePanel);
    numPlots = numPlots + 1;
    codeGenPanels{numPlots} = surfacePanel;
    comments{numPlots} = xlate( 'Plot fit with data.' );
end

if canGenerateCodeForPlot(residualsPanel)
    numPlots = numPlots + 1;
    codeGenPanels{numPlots} = residualsPanel;
    comments{numPlots} = xlate( 'Plot residuals.' );
end

if canGenerateCodeForPlot(contourPanel)
    numPlots = numPlots + 1;
    codeGenPanels{numPlots} = contourPanel;
    comments{numPlots} = xlate( 'Make contour plot.' );
end
end

function iGenerateAxesLimitCode( aFitdev, mcode )
% iGenerateDataLimitCode   Generate code that computes the limits of the axes for
% the plots
%
% This code is only required for surface, i.e., not curves, and when the
% validation data is valid.

if ~isCurve( aFitdev ) && isValidationDataValid( aFitdev )
    
    addVariable( mcode, '<xlim>', 'xlim' )
    addVariable( mcode, '<ylim>', 'ylim' )
    addBlankLine( mcode );
    addFitComment( mcode, xlate( 'Compute limits for axes.' ) );
    addFitCode( mcode, '<xlim> = [min( [<x-input>; <validation-x>] ), max( [<x-input>; <validation-x>] )];' );
    addFitCode( mcode, '<ylim> = [min( [<y-input>; <validation-y>] ), max( [<y-input>; <validation-y>] )];' );
end
end

function iGenerateContourForCurvesComment( contourPanel, mcode )
% iGenerateContourForCurvesComment will add a "NoCurveContourMessage"
% comment if the contour plot is visible, but we cannot generate code for
% the plot.
if strcmpi( contourPanel.Visible, 'on' ) && ~canGenerateCodeForPlot( contourPanel )
    addBlankLine( mcode );
    addFitComment( mcode, contourPanel.NoCurveContoursMessage.getString );
end
end

