function fns = FindFiles(globfn, varargin)

% FindFiles  Finds all files in subdirectory tree that match a wildcard input globfn
%
% fns = FindFiles(globfn, varargin)
%
% INPUTS:
%       globfn - filename to search for (you can use '*', but not '?', don't use directory names)
%       varargin PARAMETERS:
%           StartingDirectory - default '.'
%           CheckSubdirs (1/0) - default 1
% OUTPUTS:
%       fns - cell array of files found
%
% Based on matlab's dir function
% Searches all directories under the current directory
%
% ADR 1998, version L4.2, last modified '98 by ADR

% status PROMOTED
% v 4.1: added StartingDirectory parameter
% v 4.2: added CheckSubdirs parameter


%-----------------
StartingDirectory = '.';
CheckSubdirs = 1;
Extract_varargin;
pushdir(StartingDirectory);

%-----------------
thisdirfiles = dir(globfn);
fns = cell(length(thisdirfiles),1);
for iF = 1:length(thisdirfiles)
   fns{iF} = fullfile(pwd,thisdirfiles(iF).name);
end

if CheckSubdirs
   subdirs = dir;
   for iD = 1:length(dir)
      if subdirs(iD).isdir & ...
            strcmp(subdirs(iD).name,'.')==0 & ...       % is not .
            strcmp(subdirs(iD).name,'..')==0            % is not ..      
         pushdir;
         cd(subdirs(iD).name)
         % disp(['FindFiles: searching "',subdirs(iD).name,'".']);
         subfns = FindFiles(globfn);
         popdir;
         fns = [fns; subfns];
      end
   end
end

popdir;
