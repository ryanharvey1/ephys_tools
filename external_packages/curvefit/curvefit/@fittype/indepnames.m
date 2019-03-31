function names = indepnames(fun)
%INDEPNAMES Independent parameter names.
%   INDEPNAMES(FUN) returns the names of the independent parameters of the
%   FITTYPE object FUN as a cell array of strings.
%
%   See also DEPENDNAMES, FORMULA.

%   Copyright 1999-2008 The MathWorks, Inc.
%   $Revision: 1.4.2.4 $  $Date: 2009/01/23 20:37:45 $

names = cellstr(fun.indep);
