function sortedCA = sortcell(stringCA)

% sortcell  Sorts a cell array of strings
%
% sortedCA = sortcell(stringCA)
%
% INPUTS:
%       stringCA = a cell array of strings
% OUTPUTS:
%       sortedCA = sorted cell array of strings


sortedCA = cellstr(sortrows(char(stringCA)));