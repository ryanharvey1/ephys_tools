classdef ContourPlotCodeGenerator < sftoolgui.codegen.PlotCodeGenerator
    %CONTOURPLOTCODEGENERATOR   Class for generating code for contour plots
    
    %   Copyright 2009-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.4 $    $Date: 2011/03/22 18:29:55 $
    
    methods
        function obj = ContourPlotCodeGenerator()
            obj.PlotCommandGenerator = sftoolgui.codegen.SfitPlotCommandGenerator();
            obj.PlotCommandGenerator.StyleArguments = ', ''Style'', ''Contour''';
        end
    end
    
    methods(Access = protected)
        function addAxesLabels( ~, mcode )
            addFitComment( mcode, xlate( 'Label axes' ) );
            addFitCode( mcode, 'xlabel( ''<x-name>'' );' );
            addFitCode( mcode, 'ylabel( ''<y-name>'' );' );
        end
        
        function addValidationPlotCommand( obj, mcode )
            if obj.HaveLegend
                addFitCode( mcode, '<h>(end+1) = plot( <validation-x>, <validation-y>, ''bo'', ''MarkerFaceColor'', ''w'' );' );
            else
                addFitCode( mcode, 'plot( <validation-x>, <validation-y>, ''bo'', ''MarkerFaceColor'', ''w'' );' );
            end
        end
        
        function addLegendCommand( obj, mcode )
            % addLegendCommand -- Add the legend command to the generated code
            % if we have a legend
            if obj.HaveLegend
                lcg = sftoolgui.codegen.LegendCommandGenerator();

                % The first names on the legend are the fit name and the
                % fitting data name
                lcg.addName( obj.SafeFitName );
                lcg.addName( obj.SafeFittingDataName );
                
                % If there are exclusions, then they come next
                if obj.HaveExcludedData
                    lcg.addName( obj.ExcludedDataName );
                end

                % If there is validation data, then that comes last
                if obj.HaveValidation
                    lcg.addName( obj.SafeValidationDataName );
                end
                
                addCommand( lcg, mcode  );
            end
        end

    end
end
