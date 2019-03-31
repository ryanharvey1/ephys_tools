function S = LoadSpikes(tfilelist, quiet)

% LoadSpikes  Creates cell array of ts objects of spike times from list of tfiles
%
% S = LoadSpikes(tfilelist, quiet)
%
% INPUTS: 
%   tfilelist = cell array of strings, each of which is a tfile to open  
%       (incompatible with version unix3.1)
%   quiet (optional):  If quiet = 1, no messages are printed
% OUTPUTS: 
%   S = cell array such that each cell contains a ts object 
%       (timestamps which correspond to times at which the cell fired)
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status: PROMOTED


%-------------------
% Check input type
%-------------------
if nargin == 1
    quiet = 0;
end
if ~isa(tfilelist, 'cell')
   %error('LoadSpikes: tfilelist should be a cell array.');
   tfilelist = {tfilelist};
end

nFiles = length(tfilelist);

%--------------------
% Read files
%--------------------

if ~quiet
    fprintf(2, 'Reading %d files.', nFiles);
end
% for each tfile
% first read the header, the read a tfile 
% note: uses the bigendian modifier to ensure correct read format.

S = cell(nFiles, 1);
for iF = 1:nFiles
%   if ~quiet, DisplayProgress(iF, nFiles, 'Title', 'LoadSpikes'); end % REMOVED BY RYAN H 10/27/16
  tfn = tfilelist{iF};
  if ~isempty(tfn)
    tfp = fopen(tfn, 'rb','b');
    if (tfp == -1)
      warning([ 'Could not open tfile ' tfn]);
    end
    
    ReadHeader(tfp);    
    S{iF} = fread(tfp,inf,'uint32');	%read as 32 bit ints
    S{iF} = ts(S{iF});
	
    fclose(tfp);
 
  end 		% if tfn valid
end		% for all files
if ~quiet, fprintf(2,'\n'); end
