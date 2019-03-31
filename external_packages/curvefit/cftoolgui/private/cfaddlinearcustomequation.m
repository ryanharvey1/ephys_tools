function status=cfaddlinearcustomequation(arg)
% CFADDLINEARCUSTOMEQUATION Helper function for CFTOOL.

% Called by CustomLinear to create a custom equation.

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.3.2.4 $  $Date: 2008/05/01 20:12:35 $

try
  f=fittype( ...
      getValue(arg,'terms'), ...
      'coeff',getValue(arg,'coeffs'), ...
      'independent',getValue(arg,'independentVariable'), ...
      'dependent',getValue(arg,'dependentVariable'));
  opts=fitoptions(f);
catch e
  status=java.lang.String( e.message );
  return
end

% Add this fittype to the list
managecustom('set',getValue(arg,'equationName'),f,opts,getValue(arg,'equationOldName'));
status=java.lang.String('OK');

%===============================================================================
function value=getValue(arg,param)
% Get the value for the given parameter out of a list of param-value pairs
for i=1:2:length(arg)
  if isequal(arg{i},param)
    value=arg{i+1};
    return;
  end
end
