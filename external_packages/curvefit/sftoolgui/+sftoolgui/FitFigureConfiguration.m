classdef FitFigureConfiguration < sftoolgui.Configuration
    %FITFIGURECONFIG configuration file for sftool FitFigures 
    %
    %   SFTOOLGUI.FITFIGURECONFIGURATION
    %
    %   Copyright 2008 The MathWorks, Inc.
    %   $Revision: 1.1.6.2 $    $Date: 2008/10/31 05:57:47 $

    properties(SetAccess = 'private', GetAccess = 'private')
         Version = 1;
    end
    properties(SetAccess = 'public', GetAccess = 'public')
         FittingConfig;
         ResidualsConfig;
         SurfaceConfig;
         ResultsConfig;
         ContourConfig;
         Visible = 'on';
         FitUUID;
         Grid = 'on';
         LegendOn = true;
    end

    methods
          function this = FitFigureConfiguration(fID)
              this.FitUUID = fID;
              this.FittingConfig = sftoolgui.FittingConfiguration;
              this.ResidualsConfig = sftoolgui.ResidualsConfiguration;
              this.SurfaceConfig = sftoolgui.SurfaceConfiguration;
              this.ResultsConfig = sftoolgui.ResultsConfiguration;
              this.ContourConfig = sftoolgui.SurfaceConfiguration;
          end
     end
end
