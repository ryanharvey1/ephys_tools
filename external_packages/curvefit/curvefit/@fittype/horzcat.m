function c = horzcat(varargin)
%HORZCAT Horizontal concatenation of FITTYPE objects (disallowed)

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.4 $  $Date: 2010/10/08 16:37:19 $

error(message('curvefit:fittype:horzcat:catNotPermitted', class( varargin{ 1 } )));
