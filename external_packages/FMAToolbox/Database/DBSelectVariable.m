%DBSelectVariable - Select figures in current database.
%
% Select variables. Also get code used to generate them.
%
%  USAGE
%
%    [v,comments,parameters,code,date,user] = DBSelectVariable(eid,name)
%
%    eid            experiment ID (identifier string)
%    name           variable descriptive name (identifier string)
%
%  OUTPUT
%
%    v              variable
%    comments       comments
%    parameters     parameters
%    code           text of m-files used to generate the variable
%    date           date when the variable was saved
%    user           the user that saved the variable

% Copyright (C) 2007 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [v,comments,parameters,code,date,user] = DBSelectVariable(eid,name)

% Make sure MyM is installed and functional
CheckMyM;

% Check number of parameters
if nargin ~= 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help DBSelectVariable">DBSelectVariable</a>'' for details).');
end

[v,comments,parameters,c,date,user] = mym(['select v,comments,parameters,mfiles,code,date,user from variables where eid="' eid '" and name="' name '"']);
if isempty(v),
	warning(['Variable (' eid ',' name ') not found.']);
end

for i = 1:length(c{1}),
	code{i} = char(c{1}{i})';
end

v = v{1};