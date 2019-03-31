function fns = FindDirs(dirname, varargin)

% FindDirs  Finds directories in a subdirectory tree that match a wildcard input globfn
% This is usefule if you like to find e.g. all dirs called 'tfiles' in a
% directory tree. The return cell array will hold all full paths ending in
% '.../tfiles'  starting from the starting directory
%
% fns = FindFiles(globdn, varargin)
%
% INPUTS:
%       dirname - directory name to search for (EXACT match!! No wildcards like '*', etc.!!)
%       varargin PARAMETERS:
%           StartingDirectory - default '.'
%           CheckSubdirs (1/0) - default 1
% OUTPUTS:
%       fns - cell array of dirs found
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
d = dir('.');   % get all file & dirnames in current  
ix = find([d(:).isdir]);  % list of indices of all subdirs in current dir 
subdirs = d(ix(3:end));    % reduce to list of subdirs only (remove '.' and '..' dir) 
fns = {};
for iF = 1:length(subdirs)
   if strcmp(subdirs(iF).name, dirname) 
        fns{end+1,1} = [pwd filesep subdirs(iF).name];
        disp(['FindDirs found: ',fns{end,1}]);
   end
end

if CheckSubdirs
   for iD = 1:length(subdirs)
         pushdir(subdirs(iD).name);
         %disp(['FindDirs: searching "',subdirs(iD).name,'".']);
         subfns = FindDirs(dirname);
         popdir;
         fns = [fns; subfns];
   end
end

popdir;
