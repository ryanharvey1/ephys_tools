function probs = probvalues( obj )
%PROBVALUES   Problem parameter values.
%   PROBVALUES(SF) returns the values of the problem parameters of the SFIT
%   object SF as a row vector.
%
%   See also SFIT/COEFFVALUES, FITTYPE/FORMULA.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2008/10/31 05:57:03 $

if isempty( obj.fProbValues )
   probs = {};
else
   probs = [obj.fProbValues{:}];
end
