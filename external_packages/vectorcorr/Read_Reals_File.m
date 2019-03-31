function [RtrnCode, NRows, data] = Read_Reals_File(dsnIn, NSkipHdr, NCols)

%Read_Reals_File.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Open file containing ascii data to be processed
% Read entire file to count the number of rows in file
% Rewind file
% Skip past any headers in the data set
% Read the data values - must not contain any characters or strings
% Reformat array into specified number of rows and columns 
% Transpose for operation

% Input variables:
%   dsnIn: Name of input data set
%   NSkipHdr: Number of header records in dataset that are to be skipped
%   NCols: Number of columns in data array
% Output variables:
%   RtrnCode: Error return  
%       0 = OK
%       1 = Could not open or find data set
%       2 = Number of entries input does not equal product NRows*NCols 
%   NRows: Number of rows (observations) in data array 
%   data: Array that receives the input file

RtrnCode = 0;

% open file and check for existence

fid = fopen(dsnIn, 'rt');                % open file
if fid < 0
    RtrnCode = 1;
    disp('ERROR *** Dataset could not be opened or was not found')
    prt1 = ['          Dataset name specified: ',dsnIn]; disp(prt1)
    disp('JOB STOPPING'); disp(' ')
    return
end

% read entire file to find the number of records (NRows)

tf = 0;
NRows = 0;
while tf == 0
    title_skip = fgets(fid);
    NRows = NRows + 1;
    tf = feof(fid);
end

NRows = NRows - NSkipHdr;
frewind(fid);

% Start Reading

  % skip past any headers

if NSkipHdr > 0
    for iii = 1:NSkipHdr                % loop over headers 
       title_skip = fgets(fid);         % and skip  
    end
end
 
  % read records

[data, count] = fscanf(fid, '%f');      % read data into vector

fclose(fid);                            % close file

% error check

if count ~= NCols*NRows
    RtrnCode = 2;
    pt1 = sprintf(['ERROR *** Number of values input (%.0f)\n', ... 
                   '          not equal to product NRows*NCols', ...
                   ' (%.0f, %.0f)\n', ...
                   '          Likely error in NCols, or\n', ...
                   '          Characters are in data columns, or\n', ...
                   '          Incomplete record in data set\n'], ...
                   count, NRows, NCols);
    disp(pt1)
    return
end

% reformat data so that columns are variables
%   and rows are observations 

data = reshape(data, NCols, NRows);
data = transpose(data);
    
