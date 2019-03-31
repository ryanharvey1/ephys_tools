function cftoolPostEditEquationName(cftoolFits, equationName)
%CFTOOLPOSTEDITEQUATIONNAME performs update tasks after an equation name is edited
%
%   CFTOOLPOSTEDITEQUATIONNAME(CFTOOLFITS, EQUATIONNAME)
%
%   This method should be called when only the NAME of the custom
%   equation has changed. If other parts of the equation have changed
%   then you should use CFTOOLPOSTEDITEQUATION.
%
%   Also note that, unlike CFTOOLPOSTEDITEQUATION, this method will
%   perform an update to the current fit.                        FIX ME: DOES THIS WORK?
%
%   See also CFTOOLPOSTEDITEQUATION.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2007/06/14 04:55:07 $

% For each fit in the vector of fits from cftool, we need to
% -- set the new name
% -- notify cftool that fit has changed
for i = 1:length(cftoolFits)
    thisFit = cftoolFits{i};

    % We can sometimes be passed empty fits. We need to ignore these.
    if isempty( thisFit )
        continue;
    end

    % For fits to custom equations, the hint contains the name of the
    % equation.
    thisFit.hint = equationName;

    % Notify interested parties that the fit has changed
    com.mathworks.toolbox.curvefit.FitsManager.getFitsManager.fitChanged( ...
        java( thisFit ), thisFit.isGood );
end

end
