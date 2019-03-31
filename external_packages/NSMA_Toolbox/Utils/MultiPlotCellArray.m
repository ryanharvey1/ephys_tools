function MultiPlotCellArray(F, C, varargin)

% MultiPlotCellArray  Creates multiple plots by passing multiple sets of inputs to a given function
%
% MultiPlotCellArray(F, C, parameters)
%
% INPUTS:
%       F - Matlab function (string)
%       C - Cell Array, where each element = a list of inputs for function F
%       varargin PARAMETERS:
%           maxPerPage (default 64) - how many plots to put on a page
%           Name (default 'MultiPlot') - name of figure
% OUTPUTS:
%       (none)
%
% ADR 1998, version L4.2, last modified 1/18/99 by ADR

% status PROMOTED
% v 4.1 18 jan 1999 no longer uses first figures, creates own instead


%--------------------------
% Parameters

maxPerPage = 64;
Name = 'MultiPlot';
Extract_varargin;

%---------------------------
% check inputs
if ~isa(C, 'cell'); error('TYPE ERROR: C is not a cell array.'); end
if (exist(F) < 2) | (exist(F) > 6); error('TYPE ERROR: function not found.'); end

nElts = length(C);
nPages = ceil(nElts/maxPerPage);
[nL, ppL] = Subplots(min(maxPerPage, nElts));

iC = 1;
for iPage = 1:nPages
   figure('NumberTitle', 'Off', 'Name', [Name, ' -- page ', num2str(iPage)]);
   if iPage < nPages
      nThisPage = maxPerPage;
   else
      nThisPage = nElts - iC + 1;
   end  
   for iThisPage = 1:nThisPage
      subplot(nL, ppL, iThisPage);
      feval(F, C{iC});
      title(iC);
      drawnow;
      iC = iC + 1;
   end
end
