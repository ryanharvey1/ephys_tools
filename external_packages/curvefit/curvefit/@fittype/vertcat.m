function c = vertcat(varargin)
%VERTCAT Vertical concatenation of FITTYPE objects (disallowed)

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.4 $  $Date: 2010/10/08 16:37:22 $

error(message('curvefit:fittype:vertcat:catNotAllowed', class( varargin{ 1 } )));