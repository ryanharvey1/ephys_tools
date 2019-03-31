%DBInsertVariable - Save variable in current database.
%
% Save variable and code used to compute it.
%
%  USAGE
%
%    DBInsertVariable(v,eid,name,parameters,comments,mfiles)
%
%    v              variable
%    eid            experiment ID (identifier string)
%    name           descriptive name
%    comments       comment string
%    parameters     parameters
%    mfiles         cell array of m-file names (used to compute the variable)
%
%  NOTE
%
%    The pair (eid,name) must uniquely define the variable (they are used as
%    primary keys in the database).

% Copyright (C) 2007-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function DBInsertVariable(v,eid,name,parameters,comments,mfiles)

global dbUser;

% Make sure MyM is installed and functional
CheckMyM;

% Check number of parameters
if nargin < 6,
	error('Incorrect number of parameters (type ''help <a href="matlab:help DBInsertVariable">DBInsertVariable</a>'' for details).');
end

% Check parameters
if ~ischar(eid),
	error('Incorrect EID string (type ''help <a href="matlab:help DBInsertVariable">DBInsertVariable</a>'' for details).');
end
if ~ischar(name),
	error('Incorrect figure name (type ''help <a href="matlab:help DBInsertVariable">DBInsertVariable</a>'' for details).');
end
if ~ischar(comments),
	error('Incorrect comment string (type ''help <a href="matlab:help DBInsertVariable">DBInsertVariable</a>'' for details).');
end
if isempty(comments),
	comments = '-';
end
if ~ischar(parameters),
	error('Incorrect figure parameters (type ''help <a href="matlab:help DBInsertVariable">DBInsertVariable</a>'' for details).');
end
if ~iscell(mfiles),
	error('Incorrect list of m-files (type ''help <a href="matlab:help DBInsertVariable">DBInsertVariable</a>'' for details).');
end

% M-Files
for i = 1:length(mfiles),
	mfileName = which(mfiles{i});
	if isempty(mfileName),
		error(['M-File for function ''' mfiles{i} ''' not found.']);
	end
	fid = fopen(mfileName);
	code{i} = fread(fid);
	fclose(fid);
end

% Date and md5
d = datestr(now);
try
	md5 = CalcMD5(v);
catch
	warning('Could not compute MD5 (works only for numeric arrays)');
	md5 = 0;
end

% Execute SQL query command
h = mym(['insert into variables (eid,name,comments,parameters,mfile,code,date,user,md5,v) values ("{S}","{S}","{S}","{S}","{M}","{M}","{S}","{S}","{Si}","{M}")'],eid,name,comments,parameters,mfiles,code,d,dbUser,md5,v);
