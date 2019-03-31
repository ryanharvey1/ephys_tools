function c = cat(varargin)
%CAT    N-D concatenation of FITTYPE objects (disallowed)

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.4 $  $Date: 2010/10/08 16:37:15 $

error(message('curvefit:fittype:cat:catNotPermitted', class( varargin{ 1 } )))
