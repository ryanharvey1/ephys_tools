function names = dependnames(fun)
%DEPENDNAMES Dependent parameter names.
%   DEPENDNAMES(FUN) returns the names of the dependent parameters of the
%   FITTYPE object FUN as a cell array of strings.
%
%   See also FITTYPE/INDEPNAMES, FITTYPE/FORMULA.

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.2 $  $Date: 2005/03/07 17:25:59 $

names = cellstr(fun.depen);