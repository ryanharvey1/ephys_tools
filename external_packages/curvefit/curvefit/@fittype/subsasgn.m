function c = subsasgn(FITTYPE_OBJ_, varargin)
%SUBSASGN    subsasgn of fittype objects.

%   Copyright 1999-2004 The MathWorks, Inc.
%   $Revision: 1.4.2.4 $  $Date: 2010/10/08 16:37:20 $

error(message('curvefit:fittype:subsasgn:subsasgnNotAllowed', class( FITTYPE_OBJ_ )));