function str = char(obj)
%CHAR Convert FITTYPE object to character array.
%   CHAR(FUN) returns the formula for the FITTYPE object FUN.
%   This is the same as FORMULA(FUN).
%
%   See also FITTYPE, FITTYPE/FORMULA.

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.5.2.2 $  $Date: 2005/03/07 17:25:57 $

str = obj.defn;
