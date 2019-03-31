function num = numcoeffs(model)
%NUMCOEFFS Number of coefficients.
%   NUMCOEFFS(FITTYPE) returns the number of coefficients of FITTYPE.
%   
%
%   See also FITTYPE/COEFFNAMES.

%   Copyright 2001-2006 The MathWorks, Inc. 
%   $Revision: 1.2.2.2 $  $Date: 2006/12/15 19:26:17 $

num = size(model.coeff,1);
