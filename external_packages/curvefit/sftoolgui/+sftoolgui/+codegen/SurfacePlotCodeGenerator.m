classdef SurfacePlotCodeGenerator < sftoolgui.codegen.PredictablePlotCodeGenerator
    %SURFACEPLOTCODEGENERATOR   Class for generating code for surface plots
    
    %   Copyright 2009-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.3 $    $Date: 2011/03/22 18:30:05 $

    properties
        % View -- 1-by-2 vector -- view angle
        %   The angle for the view from which the observer sees the plot.
        %
        %   See also: VIEW
        View = iDefaultView;
    end

    methods
        function obj = SurfacePlotCodeGenerator()
            obj.PlotCommandGenerator = sftoolgui.codegen.SfitPlotCommandGenerator();
        end
    end
    
    methods(Access = protected)
        function str = getDefaultPredictionStyle( ~ )
            str = ', ''Style'', ''PredObs''';
        end
        
        function str = getCustomPredictionFormat( ~ )
            str = ', ''Style'', ''PredObs'', ''Level'', %s';
        end
        
        function addValidationPlotCommand( obj, mcode )
            % addValidationPlotCommand -- Add code for validation plot
            if obj.HaveLegend,
                addFitCode( mcode, '<h>(end+1) = plot3( <validation-x>, <validation-y>, <validation-z>, ''bo'', ''MarkerFaceColor'', ''w'' );' );
            else
                addFitCode( mcode, 'plot3( <validation-x>, <validation-y>, <validation-z>, ''bo'', ''MarkerFaceColor'', ''w'' );' );
            end
        end
        
        function addViewCode( obj, mcode )
            % addViewCode -- Adds code to set the view.
            %   Only display the view command if the view differs from the
            %   default 3D view.
            if any( obj.View ~= iDefaultView()  )
                codeCommand = sprintf( 'view( %.1f, %.0f );', obj.View );
                addFitCode( mcode, codeCommand );
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

function ae = iDefaultView()
% iDefaultView -- The default view angle for 3D plots
ae = [-37.5, 30];
end
