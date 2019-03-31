function thedsdb=getdsdb(varargin)
% GETDSDB is a helper function for CFTOOL

% GETDSDB Get the data set data base (singleton)

%   $Revision: 1.11.2.4 $  $Date: 2007/06/14 04:54:40 $
%   Copyright 2000-2007 The MathWorks, Inc.


thedsdb = cfgetset('thedsdb');

% Create a singleton class instance
if isempty(thedsdb)
   thedsdb = cftool.dsdb;
   cfgetset('thedsdb',thedsdb);
end


