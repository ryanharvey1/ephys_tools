function num = numargs(model)
%NUMARGS Number of argument.
%   NUMARGS(FITTYPE) returns the number of arguments of FITTYPE.
%   
%
%   See also FITTYPE/ARGNAMES.

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.2 $  $Date: 2005/03/07 17:26:16 $

num = model.numArgs;
