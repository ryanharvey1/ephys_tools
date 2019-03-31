function cella = linearexprs(model)
%LINEAREXPRS Vectorized expressions for linear coefficient matrix.
%   LINEAREXPRS(FITTYPE) returns the cell array of linear terms of FITTYPE after
%   they have been "vectorized".
%
%   See also FITTYPE/COEFFNAMES.

%   Copyright 2001-2006 The MathWorks, Inc. 
%   $Revision: 1.2.2.2 $  $Date: 2006/12/15 19:26:14 $

cella = model.Aexpr;
