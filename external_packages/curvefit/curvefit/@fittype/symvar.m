function vars = symvar(s)
%SYMVAR Determine the symbolic variables for an FITTYPE.
%   SYMVAR returns the variables for the FITTYPE object.
%
%   See also ARGNAMES.

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.2 $  $Date: 2005/03/07 17:26:45 $

vars = argnames(s);
