function header = ReadNlxHeader(fname)
% 
% Read a standard Neuralynx header and return each line as string in a cell
% array.
%
% INPUT:
%    fname ... full filename as a string
%
% OUTPUT:
%    header ... cell array of strings, each cell holding one line of the
%               header. No further parsing is done
%
% PL Jan 2005

fid = fopen(fname,'r');
if fid < 0
    [lastmsg, lasterrid] = lasterr; 
    error('Could not open file %s\n Last Error Message: %s\n Last Error ID: %d\n',fname, lastmsg,lasterrid);
end

NlxHeaderSize = 16384;   % Standard Neuralynx header size in bytes
A = fread(fid,NlxHeaderSize);
fclose(fid);

 
A = char(A(find(A)))';              % strip off NULL characters; convert double array to string
nl = strfind(A,sprintf('\n'));      % find indices of new line characters in A

header = cell(length(nl),1);
istart = 1;
for i=1:length(nl)
    iend = nl(i);
    header{i} =  deblank(A(istart:iend));
    istart = iend+1;
end
