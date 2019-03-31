classdef( Sealed ) LegendCommandGenerator < handle
    % LegendCommandGenerator   A class to generate code for legend commands
    
    %   Copyright 2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2011/03/22 18:29:59 $
    
    properties( SetAccess = private )
        % Names   Names of objects that will be displayed in the legend.
        Names = {}
    end
    
    properties( Constant, GetAccess = private )
        % Start   The first part of the legend command
        Start = 'legend( <h>';
        % Finish   The final part of the legend comamnd
        Finish = ', ''Location'', ''NorthEast'' );';
        % NameFormat  The format used to display names in the legend
        % comamnd
        NameFormat = ', ''%s''';
    end
    
    methods
        function addName( cg, name )
            % addName   Adds a name to a LegendCommandGenerator
            %
            % Note that the order that the names are added is the order
            % that they will be displayed in the legend comamnd
            cg.Names{end+1} = name;
        end
        
        function addCommand( cg, mcode )
            % addCommand   Add MATLAB code for a legend to the given 
            %   sftoolgui.codegen.MCode object.
            
            % The command format is concatenation of the start, the names
            % and the finish
            commandFormat = ['%s', repmat( cg.NameFormat, 1, length( cg.Names ) ), '%s'];
            command = sprintf( commandFormat, cg.Start, cg.Names{:}, cg.Finish );
            
            % Add the command to the mcode
            addFitCode( mcode, command );
        end
    end
end
