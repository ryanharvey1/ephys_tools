function c = isempty(obj)
%ISEMPTY   True for empty FITTYPE objects
%   ISEMPTY(F) returns true (1) if F is an empty FITTYPE object and
%   false (0) otherwise.

%   Copyright 1999-2007 The MathWorks, Inc.
%   $Revision: 1.4.2.3 $  $Date: 2007/05/23 18:37:30 $

c = logical( obj.isEmpty );
