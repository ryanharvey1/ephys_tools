function [ ft, errorStr ] = nonlinearEquationFittype( eqn, independentVars, dependentVars )
% nonlinearEquationFittype determines whether or not EQN is a valid
% equation.
%
%   EQN is string which represents the equation.
%   INDEPENDENTVARS is a cell array of strings which represent the
%   independent variables.
%   DEPENDENTVARS is a string which represents the dependent variable.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $    $Date: 2010/03/31 18:13:43 $

errorStr = '';

if ~ischar(eqn)
    errorStr = xlate('Equation must be a string.');
    ft = fittype;
    return;
end

if ~iscellstr(independentVars)
    errorStr = xlate('Independent variables must be a cell array of strings.');
    ft = fittype;
    return;
end

if ~ischar(dependentVars)
    errorStr = xlate('Dependent variable must be a string.');
    ft = fittype;
    return;
end

try
    ft = fittype( eqn, 'independent', independentVars, 'dependent', dependentVars );
catch ME
    switch ME.identifier
        case 'curvefit:fittype:missingIndVar'
            errorStr =  iHandleMissingIndVar(independentVars, ME.message);
        case 'curvefit:fittype:TooManyInputsForLibraryModel'
            errorStr = sprintf( '''%s'' is not a valid equation.', eqn);
        case 'curvefit:fittype:emptyExpression'
            errorStr = xlate('Cannot create an equation from an empty expression.' );
        otherwise
            errorStr = ME.message;
    end
    ft = fittype;
    return
end

% We do not expect calls to fittype that include independent and dependent
% inputs to return a library model that is not customnonlinear, but we will
% keep this check.
libraryModelType = type(ft);
if ~strcmp(libraryModelType, 'customnonlinear')
    errorStr = sprintf ('Library models, such as ''%s'', cannot be used as custom equations.', ...
        libraryModelType);
    ft = fittype;
end
end

function errorStr = iHandleMissingIndVar(independentVars, message)
% iHandleMissingIndVar customizes the "missingIndVar" error message. The
% original "missingIndVar" is appropriate when there is only one
% independent variable as is the case with cftool. We modify the message
% when there are 2 independent variables such as the case with sftool
% equations.

if length(independentVars) == 2
    errorStr =  sprintf( ...
        'Both independent variables %s and %s\nmust appear in the equation expression.', ...
        independentVars{1}, independentVars{2});
else
    errorStr = message;
end
end
