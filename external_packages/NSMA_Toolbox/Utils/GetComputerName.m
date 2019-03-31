function name = GetComputerName()
%
%    name = GetComputerName()
%
%    query operating system for computer name (works only on PC's with the
%    DOS environment variable COMPUTERNAME set to the computer's name
%
% PL Dec 2004

name = '';   % default return value
if strcmpi(computer,'PCWIN')
    % we are on a windows PC
    [status, result] = system('set COMPUTERNAME');
    name = deblank(result(14:end));
end