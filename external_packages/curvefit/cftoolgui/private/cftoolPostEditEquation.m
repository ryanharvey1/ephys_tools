function cftoolPostEditEquation(cftoolFits, equationName, currentFit)
% CFTOOLPOSTEDITEQUATION performs update tasks after an equation is edited
%
%   CFTOOLPOSTEDITEQUATION(CFTOOLFITS, EQUATIONNAME, CURRENTFIT)
%
%   Post Edit Equation tasks performed:
%   -- Update the name of the custom equation.
%   -- Update the fit options so that they are compatible with the new
%     equation. 
%   -- Re-calculating the fits for the new equation. Note that this not
%     done for the current fit. That computation is done from
%     CreateAFit.
%
%   CFTOOLFITS is a cell array of cftool.fit objects that need to be
%   recomputed because they all use the same custom equation and that
%   equation has changed.
%
%   EQUATIONNAME is char array with the name of the custom equation that
%   has changed.
%
%   CURRENTFIT is a cftool.fit object that represents the fit that is
%   currently selected in the table of fits in CFTOOL. This may or may
%   not be one the fits in CFTOOLFITS. This routine does NOT refit
%   CURRENTFIT. That has handled by CustomLinear (for linear custom
%   equations) or CustomGeneral (for general custom equations).
%
%   This function assumes that the only thing that has changed is the
%   equation.

% $Revision: 1.1.6.2 $  $Date: 2007/05/23 18:37:18 $
% Copyright 2006-2007 The MathWorks, Inc.

% CFTOOL may have not correctly set the global flag that indicates that
% the user has just hit the cancel button and so fitting should stop, so
% we set it now to indicate that we can continue.
cfoptstop( 0 );
% We don't need to check the global stop flag in the loop below because
% the "doFit" method will listen to it and stop fitting if it is changed
% to the stop value. Furthermore, there is a bunch of other updating
% that needs to be done and this has to happen for all fits affected by
% the equation edit even if the fitting of those is canceled. Hence we
% don't want to jump out of this loop early.

% Get the fit options for the custom equation. We will always use the
% version that has been stored with the equation rather than whatever is
% stored with the fit
equationFitOptions = managecustom( 'getopts', equationName );
    
% For each fit in the vector of fits from cftool, we need to
% -- check the options
% -- do the fit
% -- notify cftool that fit has changed
for i = 1:length(cftoolFits)
    thisFit = cftoolFits{i};

    % We can sometimes be passed empty fits. We need to ignore these.
    if isempty( thisFit )
        continue;
    end
    
    % We need to process the current fit differently for the other fits.
    thisFitIsCurrent = isequal( currentFit, thisFit );

    % Potentially the name of the custom equation has changed so we need
    % to set the hint -- for fits to custom equations, the hint contains
    % the name of the equation as it is stored in the local custom
    % equation library.
    thisFit.hint = equationName;

    % Update the fit options stored with the fit to be those stored
    % with the equation
    thisFit.fitoptions = equationFitOptions;

    %---------------------------------
    % DO THE FIT
    %---------------------------------
    % The current fit will be re-calculated via a difeerent means (see
    % CustomLinear and CustomGeneral).
    if ~thisFitIsCurrent
        thisFit.doFit( sprintf( 'custom: %s', equationName ) );

        % Notify interested parties that the fit has changed
        com.mathworks.toolbox.curvefit.FitsManager.getFitsManager.fitChanged( ...
            java( thisFit ), thisFit.isGood );
    end
end

end
