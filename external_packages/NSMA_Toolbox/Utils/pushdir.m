function DIRSTACK = pushdir(newdir)

% pushdir  Pushes current dir onto directory stack, cd's to new dir if given
%
% DIRSTACK = pushdir(newdir)
%
% INPUTS:
%       newdir - new directory that is cd'd to 
%       (if no input, stays in same dir)
% OUTPUTS:
%       DIRSTACK - global variable containing directory stack 
%                  (cell array with 1st cell = most recently added dir = top of stack)
%
% ADR 1998, version L4.1, last modified '98 by ADR

% status PROMOTED


global DIRSTACK

if isempty(DIRSTACK)
   DIRSTACK = {pwd};
else
   DIRSTACK = [{pwd} DIRSTACK];
end

if nargin == 1
cd(newdir);
end
