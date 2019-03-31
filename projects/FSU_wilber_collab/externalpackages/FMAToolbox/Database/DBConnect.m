function DBConnect(host,user,password)

%DBConnect - Connect to the database server.
%
% Open connection to MySQL database server. Connection information can
% be either passed explicitly as arguments, or read from your .my.cnf
% configuration file.
%
%  USAGE
%
%    DBConnect(host,user,password)
%
%    host               optional host
%    user               optional user
%    password           optional password
%
%  NOTE
%
%    Providing your login/password in the function call is not secure.
%    Please consider storing this sensitive information in your .my.cnf
%    file instead.
%
%  SEE
%
%    See also DBUse.
%

% Copyright (C) 2007-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global dbUser;

% Make sure MyM is installed and functional
CheckMyM;

% Try reading account info from config file
cfHost = '';
cfUser = '';
cfPassword = '';
file = fopen('~/.my.cnf');
if file ~= -1,
	while true,
		line = fgetl(file);
		if isempty(line) | line == -1, break; end
		equal = findstr('=',line);
		if isempty(equal), continue; end
		switch line(1:equal-1),
			case 'host',
				cfHost = line(equal+1:end);
			case 'user',
				cfUser = line(equal+1:end);
			case 'password',
				cfPassword = line(equal+1:end);
		end
	end
	fclose(file);
end

% Parse parameter list
if nargin < 3,
	password = cfPassword;
else
	warning('FMAToolbox:DBConnect:InsecureLogin','Providing your login/password in the function call is not secure.\nPlease consider storing this sensitive information in your .my.cnf file instead')
end
if nargin < 2,
	user = cfUser;
end
if nargin < 1,
	host = cfHost;
end

dbUser = user;

try
	h = mym('open',host,user,password);
catch e
   error(['Could not connect to MySQL. Check your login/password (' e.message ').']);
end
