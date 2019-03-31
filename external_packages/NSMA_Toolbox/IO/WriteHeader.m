function WriteHeader(fp, varargin)

% WriteHeader  Writes NSMA header
%
% WriteHeader(fp, H1, H2, H3, ...)
%
% INPUTS:
%    fp = file-pointer to file in which header is to be written
%    varargin PARAMETERS:  
%       H1, H2, H3, ... = lines (strings) to write out as header
% OUTPUTS: 
%   (none)
%
% ADR 1998, version U3.1, last modified '98 by ADR
% PL Dec 2004: corrected '\' escape charcter bug in fprintf-pwd

% status PROMOTED


% v 3.1 now accepts cell arrays as well as strings

fprintf(fp, '%%%%BEGINHEADER\n');
fprintf(fp, '%% Program: matlab\n');
fprintf(fp, [ '%% Date: ', datestr(now), '\n']);
fprintf(fp, [ '%% Directory: ', strrep(pwd,'\','\\'), '\n']);

if ~isempty(getenv('HOST'))
   fprintf(fp, [ '%% Hostname: ', getenv('HOST'), '\n']);

end

if ~isempty(getenv('USER'))
   fprintf(fp, [ '%% User: ', getenv('USER'), '\n']);

end

for iH = 1:length(varargin)

   if isa(varargin{iH}, 'cell')

      for jH = 1:length(varargin{iH})

         fprintf(fp, '%% %s\n', varargin{iH}{jH});

      end

   elseif isa(varargin{iH}, 'char')
      fprintf(fp, '%% %s\n', varargin{iH});

   else

      error('Unknown input type.');

   end
end
fprintf(fp, '%%%%ENDHEADER\n');
