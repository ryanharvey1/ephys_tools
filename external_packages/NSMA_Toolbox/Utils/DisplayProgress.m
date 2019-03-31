function DisplayProgress(iCounter, maxCounter, varargin)

% DisplayProgress  Progress bar
%
% DisplayProgress(iCounter, maxCounter, varargin)
% DisplayProgress close
%
% INPUTS:
%       iCounter = progress so far
%       maxCounter = when to stop
%       varargin PARAMETERS:
%           UseGraphics (0/1, default 1): if 1, shows bar on screen else prints to stderr
%           Title (default 'Progress so far'): title for waitbar if used
%           EOL (default 80): end of line
% OUTPUTS:
%       (none)
%
% ADR 1998, version L5.2, last modified 1/31/99 by ADR

% status PROMOTED
% v5.0 19 nov 98 Changed parameter from GraphicsTitle to Title
% v5.1 21 jan 99 if only one to counter then skip
% v5.2 31 jan 99 now allows 'close' input


% close handle
if ischar(iCounter) & strcmp(iCounter, 'close')
   global DisplayProgressHandle
   close(DisplayProgressHandle);
   clear global DisplayProgressHandle
   return
end

if maxCounter == 1; return; end

%-------------------
% parameters

UseGraphics = 1;
Title = 'Progress so far';

EOL = 80;
SoFar = 10;

Extract_varargin;

%--------------------
if UseGraphics
   global DisplayProgressHandle
   if isempty(DisplayProgressHandle)
      DisplayProgressHandle = waitbar(0, Title);
   else
      waitbar(iCounter/maxCounter);
   end
   if iCounter == maxCounter
      close(DisplayProgressHandle);
      clear global DisplayProgressHandle
   end
   drawnow;
else
   if iCounter == 1
      fprintf(2, [Title ': .']);  
   elseif iCounter == maxCounter
      fprintf(2, '\n');
   elseif rem(iCounter,EOL) == 0
      fprintf(2, '\n');
   elseif rem(iCounter,10) == 0
      fprintf(2, '!');
   else
      fprintf(2, '.');
   end
end

