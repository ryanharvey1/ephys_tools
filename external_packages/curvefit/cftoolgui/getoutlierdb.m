function theoutlierdb=getoutlierdb(varargin)
% GETOUTLIERDB is a helper function for CFTOOL

% GETOUTLIERDB Get the outlier (singleton)

%   $Revision: 1.4.2.3 $
%   Copyright 2001-2007 The MathWorks, Inc.

theoutlierdb = cfgetset('theoutlierdb');

% Create a singleton class instance
if isempty(theoutlierdb)
   theoutlierdb = cftool.outlierdb;
   cfgetset('theoutlierdb',theoutlierdb);
end


