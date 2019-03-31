function [ ft, fopts, customEquation, errorStr ] = createNewCurveFittype( fitTypeName, qualifier, optionNameValuePairs, startpoint, lower, upper)
% createNewCurveFittype creates a curve fittype and options.

%   Copyright 2010-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $    $Date: 2011/05/09 00:40:03 $
errorStr = '';
customEquation = '';

% First create the fittype
switch fitTypeName
    case 'Exponential'
        ft = iTranslateNonLinearParametricInputs('exp', qualifier);
    case 'Fourier'
        ft = iTranslateNonLinearParametricInputs('fourier', qualifier);
    case 'Gaussian'
        ft = iTranslateNonLinearParametricInputs('gauss', qualifier);
    case 'Interpolant'
        ft = iTranslateInterpolantInputs( qualifier);
    case 'Polynomial'
        ft = iTranslateNonLinearParametricInputs('poly', qualifier);
    case 'Power'
        ft = iTranslateNonLinearParametricInputs('power', qualifier);
    case 'Rational'
        ft = iTranslateNonLinearParametricInputs('rat', qualifier);
    case 'SmoothingSpline'
        ft = fittype('smoothingspline');
    case 'SumOfSine'
        ft = iTranslateNonLinearParametricInputs('sin', qualifier);
    case 'Weibull'
        ft = iTranslateNonLinearParametricInputs('weibull', qualifier);
    otherwise
        % Should be a custom equation
        customEquation = fitTypeName;
        [ft, errorStr] = iTranslateCustomInputs( fitTypeName);
end

% Now get the options
fopts = sftoolgui.util.createFitOptions(ft, optionNameValuePairs, startpoint, lower, upper);
end

function ft = iTranslateInterpolantInputs( qualifier )
% iTranslateInterpolantInputs translates curve interpolant inputs.

switch qualifier
    case 'cubic'
        ft = fittype( 'splineinterp' );
    case 'nearest'
        ft = fittype( 'nearestinterp' );
    case 'linear'
        ft = fittype( 'linearinterp' );
    case 'PCHIP'
        ft = fittype( 'pchipinterp' );
    otherwise
        warning(message('curvefit:sftoolgui:util:CreateNewCurveFittype:UnknownInterpolantFittype', qualifier));
        ft = fittype( 'linearinterp' );
end
end

function ft = iTranslateNonLinearParametricInputs(model, qualifier)
try
    ft = fittype( sprintf( '%s%s', model, qualifier ));
catch %#ok<CTCH>
    warning(message('curvefit:sftoolgui:util:CreateNewCurveFittype:UnknownNonLinearParametricFittype', model, qualifier));
    ft = fittype('poly1');
end
end

function [ft errorStr] = iTranslateCustomInputs( eqn )
[ft, errorStr] = cfswitchyard( 'nonlinearEquationFittype', eqn, {'x'},  'y' );
end
