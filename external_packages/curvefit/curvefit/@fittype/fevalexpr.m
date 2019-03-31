function args = fevalexpr(fun)
%FEVALEXPR expression to feval.
%   FEVALEXPR(FUN) returns the expression for the FITTYPE object FUN.
%
%   See also FITTYPE/ARGNAMES, FITTYPE/CHAR.

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.2 $  $Date: 2005/03/07 17:26:04 $

args = fun.expr;
