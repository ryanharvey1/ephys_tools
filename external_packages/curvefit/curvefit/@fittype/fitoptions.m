function fitoptions = fitoptions(fun)
%FITOPTIONS Get FitOptions field of the fittype.
%   FITOPTIONS(FUN) returns the FitOptions of the FITTYPE object FUN.
%
%   See also FITTYPE/FORMULA.

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.4.2.2 $  $Date: 2008/10/31 05:56:38 $

if isempty(fun.fFitoptions)
   fitoptions = {};
else
   fitoptions = copy(fun.fFitoptions);
end
