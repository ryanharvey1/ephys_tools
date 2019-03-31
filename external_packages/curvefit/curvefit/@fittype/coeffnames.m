function coeffs = coeffnames(fun)
%COEFFNAMES Coefficient names.
%   COEFFNAMES(FUN) returns the names of the coefficients of the
%   FITTYPE object FUN as a cell array of strings.
%
%   See also FITTYPE/FORMULA.

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.5.2.2 $  $Date: 2005/03/07 17:25:58 $

if isempty(fun.coeff)
   coeffs = {};
else
   coeffs = cellstr(fun.coeff);
end
