function fitdb=getfitdb(varargin)
% GETFITDB is a helper function for CFTOOL

% GETFITDB Get the fit (singleton)

% $Revision: 1.7.2.3 $  $Date: 2007/06/14 04:54:41 $
% Copyright 2000-2007 The MathWorks, Inc.


thefitdb = cfgetset('thefitdb');

% Create a singleton class instance
if isempty(thefitdb)
   thefitdb = cftool.fitdb;
end

cfgetset('thefitdb',thefitdb);
fitdb=thefitdb;
