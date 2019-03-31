function H = ReadHeader(fp)

% ReadHeader  Reads NSMA header, leaves file-read-location at end of header
%
% H = ReadHeader(fp)
%
%  INPUT: 
%      fp = file-pointer to file containing header (i.e. not filename)
%  OUTPUT: 
%      H = cell array, where each cell contains one line from the header
%
% ADR 1997, version L4.1, last modified '97 by ADR

% status: PROMOTED
% v4.1 17 nov 98 now works for files with no header


%---------------
% Get keys
beginheader = '%%BEGINHEADER';
endheader = '%%ENDHEADER';
iH = 1; H = {};
curfpos = ftell(fp);

%--------------
% go

% look for beginheader
headerLine = fgetl(fp);
if strcmp(headerLine, beginheader)
   H{1} = headerLine;
   while ~feof(fp) & ~strcmp(headerLine, endheader)     
      headerLine = fgetl(fp);
      iH = iH+1;
      H{iH} = headerLine;
   end
else % no header
   fseek(fp, curfpos, 'bof');
end
