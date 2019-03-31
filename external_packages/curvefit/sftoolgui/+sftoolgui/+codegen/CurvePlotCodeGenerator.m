classdef CurvePlotCodeGenerator < sftoolgui.codegen.PredictablePlotCodeGenerator
    % CURVEPLOTCODEGENERATOR   Class for generating code for curve plots
    
    %   Copyright 2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.1.2.1 $    $Date: 2011/07/18 00:31:02 $
    
    properties
        % HasXData -- boolean
        %   Set to true when generating code for a fit that has x-data defined.
        HasXData = true;
    end
    
    methods
        function cg = CurvePlotCodeGenerator()
            cg.PlotCommandGenerator = sftoolgui.codegen.CfitPlotCommandGenerator();
        end
    end
    
    methods( Access = protected )
        function str = getDefaultPredictionStyle( ~ )
            str = ', ''predobs''';
        end
        function str = getCustomPredictionFormat( ~ )
            str = ', ''predobs'', %s';
        end
        
        function addAxesLabels( cg, mcode )
            % addAxesLabels   Add commands for axes labels to generated code.
            %
            %   addAxesLabels( obj, mcode )
            acg = sftoolgui.codegen.CurveAxesLabelCodeGenerator();
            acg.HasXData = cg.HasXData;
            generateCode( acg, mcode );
        end
        
        function addValidationPlotCommand( cg, mcode )
            % addValidationPlotCommand -- Add code for validation plot
            if cg.HaveLegend,
                addFitCode( mcode, '<h>(end+1) = plot( <validation-x>, <validation-y>, ''bo'', ''MarkerFaceColor'', ''w'' );' );
            else
                addFitCode( mcode, 'plot( <validation-x>, <validation-y>, ''bo'', ''MarkerFaceColor'', ''w'' );' );
            end
        end
        
        function addLegendCommand( cg, mcode )
            if cg.HaveLegend
                lcg = sftoolgui.codegen.LegendCommandGenerator();
                
                % The legend command need to always start with the data
                % name
                lcg.addName( cg.SafeFittingDataName );
                
                % The next name is the optional excluded data
                if cg.HaveExcludedData
                    lcg.addName( cg.ExcludedDataName );
                end
                
                % The name of the fit always comes next
                lcg.addName( cg.SafeFitName );
                
                % If there are prediction bounds, then they are next.
                if cg.DoPredictionBounds
                    lcg.addName( sprintf( 'Lower bounds (%s)', cg.SafeFitName ) );
                    lcg.addName( sprintf( 'Upper bounds (%s)', cg.SafeFitName ) );
                end
                
                % The last name is the validation data, if there is any
                if cg.HaveValidation
                    lcg.addName( cg.SafeValidationDataName );
                end
                
                % Add the command to the mcode                
                addCommand( lcg, mcode );
            end
        end
    end
end
