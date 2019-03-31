function status=cfaddgeneralcustomequation(customgeneral,sp,lower,upper)
% CFADDGENERALCUSTOMEQUATION Helper function for CFTOOL.

% CFADDGENERALCUSTOMEQUATION is called by the custom equations 
% panel to save a custom equation.

%   Copyright 2001-2010 The MathWorks, Inc.
%   $Revision: 1.6.2.5 $  $Date: 2010/03/31 18:13:41 $

[f, status] = nonlinearEquationFittype(char(customgeneral.getEquation), ...
    {char(customgeneral.getIndependentVariable)}, ...
    char(customgeneral.getDependentVariable));

if isempty(status)
    opts=fitoptions(f);
    opts.StartPoint=sp;
    opts.Lower=lower;
    opts.Upper=upper;
    % Add this fittype to the list
    managecustom('set',char(customgeneral.getEquationName),f,opts,char(customgeneral.getEquationOldName));
    status=java.lang.String('OK');
end
