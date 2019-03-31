classdef AbstractPlotCommandGenerator < handle
    %AbstractPlotCommandGenerator   Abstract class for generating code to call the
    % plot method of a fit object
    
    %   Copyright 2011-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2011/03/22 18:29:53 $
    
    properties
        HaveValidation = false;
        HaveLHS = false;
        HaveExcludedData = false;
        % StyleArguments -- string
        %
        %   A string of parameter-value pairs that describe the plot style.
        %   This property should set this to an appropriate value to get
        %   the plot that is needed.
        StyleArguments = '';
    end
    
    methods( Abstract )
        addPlotCommand( cg, mcode )
    end
    
end



