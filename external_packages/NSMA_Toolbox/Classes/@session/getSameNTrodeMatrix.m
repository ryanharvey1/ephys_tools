function sameNtMatrix = getSameNTrodeMatrix(ses)
% sameNtMatrix = getSameNTrodeMatrix(ses)
% 
% for a given session object, returns a sameNtMatrix
%
%  The same-NTrode-Matrix is a N*N matrix with ones for cellpairs from different probes and
%  zeros for cellpairs from same probe from a given list of tfiles with naming
%  convention secified in ses.sameTTfmt. The ordering of matrix is the same as the ordering returned 
%  by the session methods getTFileList(ses,epoch), getSpikeMatrix(ses,epoch) and getNumCells(ses).
%
%  ses.sameTTfmt needs to be a string of the form 'fixedNamePart%##' where the '%##' are placeholders for
%  a letter ('%') and digits ('##') used to build up an electrode or tetrode identifier string in t-file names. 
%  The 'fixedNamePart' is any string before the N-trode number, fixed for all t-files in the session.
%  The '%' character denotes the position of a single letter (usually A-L) describing the row number of an 
%  electrode in Warp-type drives. For hyperdrives this letter character placeholder is usually not part of a naming
%  convention an can be ommitted. 
%  The '##' part denotes the location of the n-trode number (which can be a single or several digits)
%  in the t-filename. After the N-trode number can be any sequence of characters where the first
%  character after the tetrode number must not be a digit (otherwise it would be interpreted as part of 
%  the n-trode number).
%  
%
% PL April 2004

% get tfile names
tfiles = getTFileList(ses,'any');
nC = length(tfiles);

% get position of n-trode letter
lett_indx = findstr(ses.sameTTfmt,'%');
if length(lett_indx) > 1
    error('getSameNTrodeMatrix requires ses.sameTTfmt to have 0 or 1 characters "%".');
end%if

% get position of n-trode number
tt_indx = findstr(ses.sameTTfmt,'##');
if length(tt_indx) ~= 1
    error('getSameNTrodeMatrix requires ses.sameTTfmt needs to be a string of the form "fixedNamePart%##".');
end%if

% extract nTrode number from each filename and store in nTrodeNum list parralel to tfiles list
nTrodeNum = cell(nC,1);
for iC = 1:nC
    [dn, fn, ext] = fileparts(tfiles{iC});
    
    % extract n-trode letter identifier
    if isempty(lett_indx)
        nTrodeLett = '';
    else
        [dn, fn, ext] = fileparts(tfiles{iC});
        nTrodeLett = fn(lett_indx);
    end    
    % count number of digits in n-trode number until first non-digit character is encountered
    ndigits = 0;
    while 1
        ix = tt_indx + ndigits;
        char = fn(ix);
        if isempty(findstr(char,'0123456789'))
            break;
        else
            ndigits = ndigits+1;
        end
        if (tt_indx+ndigits) > length(fn)
           break;
        end
    end
    nTrodeNum{iC} = [nTrodeLett fn(tt_indx:tt_indx+ndigits-1)];
end


% construct sameNtMatrix
sameNtMatrix = ones(nC,nC);
for iC = 1:nC
    for jC = iC:nC
        if length(nTrodeNum{iC}) ~= length(nTrodeNum{jC})
            sameNtMatrix(iC,jC) = 1;
            sameNtMatrix(jC,iC) = sameNtMatrix(iC,jC);   
        else
            sameNtMatrix(iC,jC) = ~strcmpi(nTrodeNum{iC},nTrodeNum{jC});
            sameNtMatrix(jC,iC) = sameNtMatrix(iC,jC);   
        end
    end
end

