function [eeg, sFreq] = ReadCR_tsd(fname, start_ts, end_ts, output_sFreq)

% ReadCR_tsd  Reads a CSC file (NT format) and returns a tsd
%
% [eeg, sFreq] = ReadCR_tsd(fname, start_ts, end_ts, output_sFreq)
%
% INPUTS:
%       fname = full filename of Cheetah_NT CSC*.dat file 
%       start_ts, end_ts (optional) = start and end timestamps of portion to read
%       output_sFreq (optional) = frequency at which to resample data
% OUTPUTS:
%       eeg = tsd of the csc data
%       sFreq = sampling frequency of data in eeg output
%
% ADR 1999, version 1.0, last modified '01 by SC

% status PROMOTED
% cowen Sat Jul  3 14:59:47 1999
% lipa  modified for NT   Jul 18 1999
% o Got rid of the diplay progress and rounded the timestamps
% o Fixed the dT to be in timestamps 
% o Added a dummy value to the end of cr.ts
% o Made the for look go to nBlocks and not nBlocks -1 
% cowen 2001
% o returns sampling freq.
%       ReadCR_nt returns 2 arrays and 1 double of the form...
%       ts = nrec x 1 array of the timestamps that start each block.
%       cr = nrec x 512 array of the data
%       sFreq = sampling frequency of the data in eeg
% cowen modified to use the partial load version of ReadCR_nt


if nargin >= 3    
    [ts,cr,sFreq] = ReadCR(fname,start_ts, end_ts);  %  timestamps ts are in 0.1 milliseconds units!!!!!
elseif nargin == 1
    [ts,cr,sFreq] = ReadCR(fname);  %  timestamps ts are in 0.1 milliseconds units!!!!!
elseif nargin == 2
    error('Invalid number of inputs')
end

if nargin < 4
    output_sFreq = [];
end

nBlocks = size(cr,1);
% Reuse cr to save space.
cr=reshape(cr',1,length(cr(:)));
blockSize = 512;
dT = 10000/sFreq; % in tstamps
TIME = zeros(size(cr));
ts = [ts;ts(end) + 512*dT];     
for iBlock = 1:(nBlocks)
  %DisplayProgress(iBlock, nBlocks-1);
  TIME((blockSize * (iBlock-1) + 1):(blockSize * iBlock)) = ...
      linspace(ts(iBlock), ts(iBlock+1) - dT, blockSize);
end
% Create the tsd object. Perhaps a waste.
% pack % this command is deprecated in R18

if isempty(output_sFreq)
    eeg = tsd(TIME', cr');
else
    st = TIME(1);
    et = TIME(end);
    n_out = round(((et-st)/10000)*output_sFreq);
    sFreq = output_sFreq;
    eeg  = tsd(linspace(st,et,n_out)', interp1(TIME, cr, linspace(st,et,n_out))');
end
