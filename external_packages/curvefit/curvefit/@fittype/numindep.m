function n = numindep(fun)
%NUMINDEP Number of independent parameter names.
%   NUMINDEP(FUN) returns the number of independent parameters of the FITTYPE
%   object FUN.
%
%   See also INDEPNAMES.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2008/10/31 05:56:44 $

n = size( fun.indep, 1 );
