classdef PlotCodeGenerator < handle
    %PLOTCODEGENERATOR   Base class for plot code generators
    
    %   Copyright 2009-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.5 $    $Date: 2011/03/22 18:30:01 $
    
    properties
        % FitName -- string
        %   Name of the fit to be used in the legend
        FitName = 'My Fit 1';
        
        % FittingDataName -- string
        %   Name of the fitting data to be used in the legend
        FittingDataName = 'z vs. x, y';
        
        % ValidationDataName -- string
        %   Name of the validation data to be used in the legend
        ValidationDataName = 'zv vs. xv, yv';
        
        % HaveValidation -- boolean
        %   Set to true when there is validation data that should be plotted.
        HaveValidation = false;
        
        % GridState -- 'on' or 'off'
        %   Indicates that the grid should or should not be plotted
        GridState = 'on';
        
        % HaveLegend -- boolean
        %   Set to true when code for a legend should be generated.
        HaveLegend = false;
        
        % HaveExcludedData -- boolean
        %   Set to true when there is excluded data
        HaveExcludedData = false;
    end
    properties(SetAccess = protected, GetAccess = protected)
        % PlotCommandGenerator -- AbstractPlotCommandGenerator
        PlotCommandGenerator = [];
    end
    properties(SetAccess = private, GetAccess = protected, Dependent)
        % SafeFitName -- string
        %   A "safe" version of obj.FitName. This is safe in the sense that in
        %   can be inserted into code with out causing syntax errors. For
        %   example, any quotes (') in obj.FitName will be replaced bu double
        %   quotes ('').
        SafeFitName
        
        % SafeFittingDataName -- string
        %   A "safe" version of obj.FittingDataName.
        SafeFittingDataName
        
        % SafeFittingDataName -- string
        %   A "safe" version of obj.ValidationDataName.
        SafeValidationDataName
        
        % ExcludedDataName -- string
        %   Name of any excluded data
        ExcludedDataName
    end
    
    methods
        function name = get.SafeFitName( obj )
            name = sftoolgui.codegen.stringLiteral( obj.FitName );
        end
        function name = get.SafeFittingDataName( obj )
            name = sftoolgui.codegen.stringLiteral( obj.FittingDataName );
        end
        function name = get.SafeValidationDataName( obj )
            name = sftoolgui.codegen.stringLiteral( obj.ValidationDataName );
        end
        function name = get.ExcludedDataName( obj )
            name = sprintf( 'Excluded %s', obj.SafeFittingDataName );
        end
        
        function obj = PlotCodeGenerator()
            obj.PlotCommandGenerator = sftoolgui.codegen.SfitPlotCommandGenerator();
        end
    end

    methods( Sealed )
        function generateMCode( obj, mcode )
            % generateMCode -- Generate MATLAB code and add it to the given
            %   sftoolgui.codegen.MCode object.
            
            % Register a handle for the legend?
            addLegendVariable( obj, mcode )
            % Plot command
            addPlotCommand( obj, mcode );
            % Plot Validation Data?
            addValidationPlotBlock( obj, mcode )
            % Add a Legend?
            addLegendCommand( obj, mcode )
            % Axes labels
            addAxesLabels( obj, mcode );
            % Grid?
            addGridCommand( obj, mcode );
            % View
            addViewCode( obj, mcode );
        end
    end

    methods(Access = private)
        function addLegendVariable( obj, mcode )
            if obj.HaveLegend
                addVariable( mcode, '<h>', 'h' );
            end
        end
        
        function addPlotCommand( obj, mcode )
            % addPlotCommand -- Add the main plot command to the generated code.
            pcg = obj.PlotCommandGenerator;
            
            % If we have a legend then we need a LHS from the plot command
            pcg.HaveLHS = obj.HaveLegend;
            
            pcg.HaveValidation = obj.HaveValidation;
            pcg.HaveExcludedData = obj.HaveExcludedData;
            
            addPlotCommand( pcg, mcode );
        end
        
        function addValidationPlotBlock( obj, mcode )
            % addValidationPlotBlock -- Add a block of code that plots
            % validation data.
            %
            %   To change the style of the plot for a subclass, overload the
            %   "addValidationPlotCommand" method.
            %
            %   See also addValidationPlotCommand.
            if obj.HaveValidation
                addFitComment( mcode, xlate( 'Add validation data to plot.' ) );
                addFitCode( mcode, 'hold on' );
                addValidationPlotCommand( obj, mcode );
                addFitCode( mcode, 'hold off' );
            end
        end
                
        function addGridCommand( obj, mcode )
            % addGridCommand -- Add the GRID command to the generated code
            gridCommand = sprintf( 'grid %s', obj.GridState );
            addFitCode( mcode, gridCommand );
        end
    end
    
    methods(Abstract, Access = protected)
        % addValidationPlotCommand -- Overload this method to change the
        % plot command for validation data.
        addValidationPlotCommand( obj, mcode )
        
        % addLegendCommand -- Add the legend command to the generated code
        % if we have a legend
        addLegendCommand( obj, mcode )
    end
    methods(Access = protected)
        function addAxesLabels( ~, mcode )
            % addAxesLabels -- Overload this method to set the axes labels.
            %   By Default the labels of the x-, y- and z-axes are set to be the
            %   same name as the corresponding fitting variable.
            %
            %   >> addAxesLabels( obj, mcode )
            addFitComment( mcode, xlate( 'Label axes' ) );
            addFitCode( mcode, 'xlabel( ''<x-name>'' );' );
            addFitCode( mcode, 'ylabel( ''<y-name>'' );' );
            addFitCode( mcode, 'zlabel( ''<z-name>'' );' );
        end
        
        function addViewCode( ~, ~ )
            % addViewCode -- Overload this method to add code to set the view.
            %   By Default, no view code is added
            %
            %   >> addViewCode( obj, mcode );
        end
    end
end
