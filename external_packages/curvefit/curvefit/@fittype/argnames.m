function args = argnames(fun)
%ARGNAMES Argument names.
%   ARGNAMES(FUN) returns the names of the input arguments of the
%   FITTYPE object FUN as a cell array of strings.
%
%   See also FITTYPE/FORMULA.

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.2 $  $Date: 2005/03/07 17:25:54 $

args = cellstr(fun.args);
