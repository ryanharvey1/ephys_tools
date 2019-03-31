%DBSelectVariables - Select figures in current database.
%
% Select variables. Also get code used to generate them.
%
%  USAGE
%
%    [v,comments,parameters,code,date,user] = DBSelectVariables(query)
%
%    query          optional list query (WHERE clause; see Example)
%
%  OUTPUT
%
%    v              variables
%    comments       comments
%    parameters     parameters
%    code           text of m-files used to generate the variables
%    date           date when the variables were saved
%    user           the user that saved the variables
%
%  EXAMPLE
%
%    Get all variables from experiment "experiment1", the name of which starts
%    with "raster" (for details, see an SQL manual):
%
%    [v,comments,parameters,code,date,user] = ...
%      DBSelectVariables('eid="experiment1" and name like "raster%"');

% Copyright (C) 2007-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [v,comments,parameters,code,date,user] = DBSelectVariables(query)

% Make sure MyM is installed and functional
CheckMyM;

% Edit query
if nargin == 0,
	query = '';
end
query = strtrim(query);
query = regexprep(query,'^where','');
if ~isempty(query), query = [' where ' query]; end

% Query database
[v,comments,parameters,c,date,user] = mym(['select v,comments,parameters,mfiles,code,date,user from variables' query]);
if isempty(v),
	warning(['No variables match (' query ').']);
end

for i = 1:length(c),
	for j = 1:length(c{i}),
		code{i}{j} = char(c{i}{j})';
	end
end

