classdef CurveAxesLabelCodeGenerator < handle
    % CurveAxesLabelCodeGenerator   Class for generating code for axes labels for
    % plots of curve fits.
    
    %   Copyright 2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.1 $    $Date: 2011/07/18 00:31:01 $
    
    properties
        % HasXData -- boolean
        %   Set to true when generating code for a fit that has x-data defined.
        HasXData = true;
    end
    
    methods( Access = public )
        function generateCode( cg, mcode )
            % generateCode   Generate code for axes labels for plots of curve fits.
            %
            %   generateCode( obj, mcode )
            addFitComment( mcode, xlate( 'Label axes' ) );
            if cg.HasXData
                addFitCode( mcode, 'xlabel( ''<x-name>'' );' );
            end
            addFitCode( mcode, 'ylabel( ''<y-name>'' );' );
        end
    end
end
