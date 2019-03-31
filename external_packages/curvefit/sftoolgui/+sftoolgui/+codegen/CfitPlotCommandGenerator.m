
%   Copyright 2011 The MathWorks, Inc.

classdef( Sealed ) CfitPlotCommandGenerator < sftoolgui.codegen.AbstractPlotCommandGenerator
    % CfitPlotCommandGenerator   A class for generating code to call the
    % plot method of CFIT
    
    methods
        function addPlotCommand( cg, mcode )
            % addPlotCommand   Add the command for plotting CFIT objects to
            % the generated code.
            
            plotCommandStart = 'plot( <fo>, <x-input>, <y-input>';
            plotCommandEnd = ' );';
            
            addFitCode( mcode, [
                plotCommandLHS( cg ), ...
                plotCommandStart, ...
                exclusionPlotOption( cg ), ...
                cg.StyleArguments, ...
                plotCommandEnd
                ] );
        end
    end
    
    methods( Access = private )
        function str = plotCommandLHS( cg )
            if cg.HaveLHS
                str = '<h> = ';
            else
                str = '';
            end
        end
        
        function str = exclusionPlotOption( cg )
            if cg.HaveExcludedData
                str = ', <ex>';
            else
                str = '';
            end
        end
    end
end
