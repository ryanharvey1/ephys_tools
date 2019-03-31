
%   Copyright 2011 The MathWorks, Inc.

classdef PredictablePlotCodeGenerator < sftoolgui.codegen.PlotCodeGenerator
    % PredictablePlotCodeGenerator   Abstract class for PlotCodeGenerator
    % classes that support display of prediction bounds.
    
    properties(Dependent)
        % DoPredictionBounds -- boolean
        %   Set to true to include code for plotting prediction bounds
        DoPredictionBounds;
        
        % PredictionLevel -- scalar \in (0, 1.0)
        %   The prediction level for the surface
        PredictionLevel
    end
    
    properties(Access = protected)
        % PDoPredictionBounds -- private storage of DoPredictionBounds
        PDoPredictionBounds = false;
        
        % PPredictionLevel -- private storage of PredictionLevel
        PPredictionLevel = iDefaultPredictionLevel()
    end
    
    methods
        function obj = set.DoPredictionBounds( obj, tf )
            obj.PDoPredictionBounds = tf;
            obj = updatePlotStyleArgs( obj );
        end
        function tf = get.DoPredictionBounds( obj )
            tf = obj.PDoPredictionBounds;
        end
        
        function obj = set.PredictionLevel( obj, level )
            obj.PPredictionLevel = level;
            obj = updatePlotStyleArgs( obj );
        end
        function level = get.PredictionLevel( obj )
            level = obj.PPredictionLevel;
        end
    end
    
    methods(Access = private)
        function obj = updatePlotStyleArgs( obj )
            % If we are doing prediction bounds then we need to set the style
            % appropriately.
            if obj.PDoPredictionBounds
                if obj.PPredictionLevel == iDefaultPredictionLevel()
                    styleArguments = getDefaultPredictionStyle( obj );
                else
                    styleArguments = sprintf( getCustomPredictionFormat( obj ), ...
                        num2str( obj.PredictionLevel ) );
                end
            else
                styleArguments = '';
            end
            
            obj.PlotCommandGenerator.StyleArguments = styleArguments;
        end
    end
    
    methods( Abstract, Access = protected )
        str = getDefaultPredictionStyle( obj )
        
        % getCustomPredictionFormat
        %
        %   Should contain a %s indicating where the level should be
        %   inserted.
        str = getCustomPredictionFormat( obj )
    end
end

function level = iDefaultPredictionLevel()
% iDefaultPredictionLevel -- The default value for the prediction interval level
defaultPredictionOpts = curvefit.PredictionIntervalOptions;
level = defaultPredictionOpts.Level;
end
