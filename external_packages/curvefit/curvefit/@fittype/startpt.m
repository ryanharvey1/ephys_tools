function num = startpt(model)
%STARTPT function to compute start point.
%   STARTPT(FITTYPE) returns the function handle of a function
%   to compute a starting point for FITTYPE based on xdata and ydata.
%   
%
%   See also FITTYPE/FORMULA.

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.3.2.2 $  $Date: 2008/10/31 05:56:47 $

num = model.fStartpt;
