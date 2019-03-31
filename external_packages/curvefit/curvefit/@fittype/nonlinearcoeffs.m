function cella = nonlinearcoeffs(model)
%NONLINEARCOEFFS array of indices of nonlinear coefficients.
%   NONLINEARCOEFFS(FITTYPE) returns the array of indices of 
%   nonlinear coefficients of FITTYPE.
%
%   See also FITTYPE/COEFFNAMES.

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.2.2.3 $  $Date: 2008/10/31 05:56:43 $

cella = model.fNonlinearcoeffs;
