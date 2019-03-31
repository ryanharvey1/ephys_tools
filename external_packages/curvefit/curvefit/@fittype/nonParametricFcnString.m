function line = nonParametricFcnString(obj, variable, arglist, coefficient)
%NONPARAMETRICFCNSTRING Make the string that describes a non-parametric function
%
%   LINE = NONPARAMETRICFCNSTRING(OBJ, VARIABLE, ARGLIST, COEFFICIENT)

%   Copyright 2008-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/02/17 18:45:45 $ 

lhs = leftHandSideForFcnString( variable, arglist );
definition = definitionFromLibname( numindep( obj ), obj.fType );
line = sprintf( definition, lhs, coefficient );
end

