function F = ReadFileList(fn)

% ReadFileList  Reads list of files from ascii file, creates cell array of filenames
%
% F = ReadFileList(fn)
%
% INPUTS: 
%       fn - an ascii file of filenames, 1 filename per line
% OUTPUTS:
%       F - a cell array of filenames suitable for use in programs
%           such as LoadSpikes
%
%
% ADR 1998, version L4.1, last modified '98 by ADR

% status: PROMOTED
% v4.1 added ReadHeader, now can handle files with headers


[fp,errmsg] = fopen(fn, 'rt');
if (fp == -1)
   error(['Could not open "', fn, '". ', errmsg]);
end

ReadHeader(fp);
ifp = 1;
while (~feof(fp))
   F{ifp} = fgetl(fp);
   ifp = ifp+1;
end
fclose(fp);

F = F';

