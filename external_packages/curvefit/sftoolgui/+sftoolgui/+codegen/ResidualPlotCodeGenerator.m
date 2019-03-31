classdef ResidualPlotCodeGenerator < sftoolgui.codegen.PlotCodeGenerator
    %RESIDUALPLOTCODEGENERATOR   Class for generating code for residual plots
    
    %   Copyright 2009-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.5 $    $Date: 2011/03/22 18:30:03 $
    
    methods(Sealed, Access = protected)
        function addValidationPlotCommand( obj, mcode )
            
            % If we have a legend, then we need to catch the handle to the
            % validation line in the LHS (left hand side) of the command.
            if obj.HaveLegend
                lhs = '<h>(end+1) = ';
            else
                lhs = '';
            end
            
            % The fit code is just the concatenation of the LHS and the RHS
            addFitCode( mcode, [lhs, getValidationPlotRHS( obj )] );
        end
    end
    
    methods( Abstract, Access = protected )
        % getValidationPlotRHS   The RHS of the command to plot validation
        % data
        %
        % Sub-classes need to overload this method and use it to the define
        % the RHS (right hand side) of the command used to plot validation
        % data.
        str = getValidationPlotRHS( obj )
    end
end
