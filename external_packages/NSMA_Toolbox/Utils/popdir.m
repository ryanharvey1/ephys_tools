function NEWDIR = popdir(all);

% popdir  Pops top directory off current directory stack, cd's to it
%
% NEWDIR = popdir
% NEWDIR = popdir('all')
%
% INPUTS:
%       'all' - optional: cd's to bottom dir in stack, empties stack
% OUTPUTS:
%       NEWDIR - dir that is cd'd to
%
% directory stack maintained in global variable 
%   DIRSTACK (cell array with 1st cell = most recently added dir = top of stack)
%
% ADR 1998, version L4.1, last modified '98 by ADR

% status PROMOTED
% allows popdir all now


global DIRSTACK

switch nargin
case 1
   if strcmp(all, 'all')
      if ~isempty(DIRSTACK)
         NEWDIR = DIRSTACK{length(DIRSTACK)};
         cd(NEWDIR);
         DIRSTACK = {};
      end
   else 
      error('Unknown input arguments.');
   end
case 0
   if isempty(DIRSTACK)
      disp('Directory stack empty.');
   else
      NEWDIR = DIRSTACK{1};
      cd(NEWDIR);
      DIRSTACK = DIRSTACK(2:length(DIRSTACK));
   end 
otherwise
   error('Unknown input arguments.');
end
