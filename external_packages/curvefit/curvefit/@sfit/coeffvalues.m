function coeffs = coeffvalues( obj )
%COEFFVALUES   Coefficient values.
%   COEFFVALUES(SF) returns the values of the coefficients of the
%   SFIT object SF as a row vector.
%
%   See also SFIT/PROBVALUES, FITTYPE/FORMULA.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2008/10/31 05:56:56 $

if isempty( obj.fCoeffValues )
   coeffs = {};
else
   coeffs = [obj.fCoeffValues{:}];
end
