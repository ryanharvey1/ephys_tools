function cella = linearterms(model)
%LINEARTERMS Cell array of linear terms to form linear coefficient matrix.
%   LINEARTERMS(FITTYPE) returns the cellarray of linear terms of FITTYPE.
%
%   See also FITTYPE/COEFFNAMES.

%   Copyright 2001-2006 The MathWorks, Inc. 
%   $Revision: 1.3.2.2 $  $Date: 2006/12/15 19:26:15 $

cella = model.Adefn;
