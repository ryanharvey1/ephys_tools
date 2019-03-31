function cfsaveanalysis(varnames, data)
%CFSAVEANALYSIS saves analysis information to the workspace
%This is a helper function for the Curve Fitting tool

%   Copyright 2001-2005 The MathWorks, Inc.
%   $Revision: 1.6.2.2 $  $Date: 2005/03/07 17:25:30 $

results = cfgetset('analysisresults');
if isempty(results) || ~isstruct(results) || ~isfield(results,'xi') ...
                    || isempty(results.xi)
   uiwait(warndlg('No analysis results available.',...
                  'Save to Workspace','modal'));
else
   assignin('base',varnames{1},results);
end


