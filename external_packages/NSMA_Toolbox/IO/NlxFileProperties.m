function props = NlxFileProperties(fname,ftype)
%
%  Get File properties of a Neuralynx data file and return a properties
%  struct props.
%
% INPUT: 
%   fname   .... full filename of Neuralynx file
%   ftype   .... 'SE' = Nlx single electrode file (Windows Cheetah file) 
%                'ST' = Nlx stereotrode file (Windows Cheetah file) 
%                'TT' = Nlx tetrode file (Windows Cheetah file) 
%                'CSC' = Nlx Continous Sample file (Windows Cheetah file) 
%                'EV' = Nlx event file (Windows Cheetah file) 
%                'VT' = Nlx video file (Windows Cheetah file) 
% OUTPUT:
%    props struct with fields:
%     props.fileType =  same as input ftype
%     props.recordLength_Bytes =  number of bytes in each record;
%     props.header = Cell array of strings; one string per non-empty line in the header;
%     props.headerSize_Bytes =  Standard Neuralynx header size in bytes
%     props.numRecords = number of records in file;
%     props.firstLastTS_usec = 2*1 array with first and last timestamps in native Neuralynx (1 microSec) untis;
%
%
% Warning: Work in progress .... 
% PL Jan 2005

if nargin == 1
    % infer filetype from filename extension
    [dd,fn,ext] = fileparts(fname);
    switch lower(ext)
        case '.nse'
            ftype = 'SE';
        case '.nst'
            ftype = 'ST';
        case '.ntt'
            ftype = 'TT';
        case '.ncs'
            ftype = 'CSC';
        case '.nev'
            ftype = 'EV';
        case '.nvt'
            ftype = 'VT';
        case '.dat'
            error('Neuralynx .dat files require that you specify the filetype explicitly in second argument. Read help NlxFileProperties!');
        otherwise
            error('Unknown Neuralynx file type')
    end
        
end


NlxHeaderSize = 16384;   % Standard Neuralynx header size in bytes

props.fileType = '';
props.recordLength_Bytes = 0;
props.header = ReadNlxHeader(fname);
props.headerSize_Bytes = NlxHeaderSize;   % Standard Neuralynx header size in bytes
props.numRecords = 0;
props.firstLastTS_usec = [];

% open file and position filepointer after header:
fid = fopen(fname,'r');
if fid < 0
    [lastmsg, lasterrid] = lasterr; 
    error('Could not open file %s\n Last Error Message: %s\n Last Error ID: %d\n',fname, lastmsg,lasterrid);
end
fseek(fid,0,'eof');    % goto end of file
fileSizeBytes = ftell(fid);
    if fileSizeBytes == -1
        error('Could not determine the size of file %s',fname);
    end
fseek(fid,NlxHeaderSize,'bof');   % goto end of header


% get numRecords and firstLastTS_usec
switch ftype
    
    case 'SE'
        props.fileType = 'SE';
        recSize = 112;  % record size in bytes
        tsOffset = 0;
        props.recordLength_Bytes = recSize;
        [numRecs, t0, t1] = GetNumRecsT0T1(fid,recSize,tsOffset,fileSizeBytes,NlxHeaderSize);
        props.numRecords = numRecs;
        props.firstLastTS_usec = [t0 t1];

    case 'ST'
        props.fileType = 'ST';
        recSize = 8+10*4+32*2*2;  % record size in bytes
        tsOffset = 0;
        props.recordLength_Bytes = recSize;
        [numRecs, t0, t1] = GetNumRecsT0T1(fid,recSize,tsOffset,fileSizeBytes,NlxHeaderSize);
        props.numRecords = numRecs;
        props.firstLastTS_usec = [t0 t1];

    case 'TT'
        props.fileType = 'TT';
        recSize = 8+10*4+32*4*2;  % record size in bytes
        tsOffset = 0;
        props.recordLength_Bytes = recSize;
        [numRecs, t0, t1] = GetNumRecsT0T1(fid,recSize,tsOffset,fileSizeBytes,NlxHeaderSize);
        props.numRecords = numRecs;
        props.firstLastTS_usec = [t0 t1];

    case 'CSC'
        props.fileType = 'CSC';
        recSize = 8+3*4+512*2;  % record size in bytes
        tsOffset = 0;
        props.recordLength_Bytes = recSize;
        [numRecs, t0, t1] = GetNumRecsT0T1(fid,recSize,tsOffset,fileSizeBytes,NlxHeaderSize);
        props.numRecords = numRecs;
        props.firstLastTS_usec = [t0 t1];

    case 'EV'
        props.fileType = 'EV';
        recSize = 3*2+8+3*2+9*4+128*1;  % record size in bytes
        tsOffset = 3*2;
        props.recordLength_Bytes = recSize;
        [numRecs, t0, t1] = GetNumRecsT0T1(fid,recSize,tsOffset,fileSizeBytes,NlxHeaderSize);
        props.numRecords = numRecs;
        props.firstLastTS_usec = [t0 t1];

    case 'VT'
        props.fileType = 'VT';
        recSize = 1828;  % record size in bytes
        tsOffset = 3*2;
        props.recordLength_Bytes = recSize;
        [numRecs, t0, t1] = GetNumRecsT0T1(fid,recSize,tsOffset,fileSizeBytes,NlxHeaderSize);
        props.numRecords = numRecs;
        props.firstLastTS_usec = [t0 t1];

        
    otherwise
        error('ftype = %s is not implemented!')
end


fclose(fid);
return


%%-------------------------------------------------------------------------
function [numRecs, t0, t1] = GetNumRecsT0T1(fid,recSize,tsOffset,fileSizeBytes,NlxHeaderSize)
% 
numRecs = (fileSizeBytes - NlxHeaderSize)/recSize;
junkBytes = 0;
if (numRecs - floor(numRecs)) > 0
    junkBytes = (fileSizeBytes - NlxHeaderSize) - floor(numRecs)*recSize;
    warning('Number of Records %f is not an integral number. Ignoring last %d bytes in file!',junkBytes);
end
numRecs = floor(numRecs);
% read first and last timestamp
fseek(fid,NlxHeaderSize+tsOffset,'bof');   % goto end of header
t0 = fread(fid,1,'int64');   % fist timestamp
fseek(fid,-recSize-junkBytes+tsOffset,'eof');   % goto last record
t1 = fread(fid,1,'int64');   % last timestamp
return
