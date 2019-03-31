function [avgEEG, errEEG, EEGts, EEGsamples] = PETH_EEG(T, eeg_fname, tlow_msec, thigh_msec)

% PETH_EEG  Calculates PETH of EEG (rather than spike) data
%
% [avgEEG, errEEG, ts, EEGsamples] = PETH_EEG(T, eeg_fname, tlow_msec, thigh_msec)  
%
% INPUTS: 
%   T - array of timestamps (0.1msec units) of "centers" of the PETH
%   eeg_fname - path/filname of EEG file
%   tlow_msec - lower edge T-tlow_msec of PETH window ( in msec  units)
%   thigh_msec - upper edge T+thigh_msec of PETH window ( in msec  units)
% OUTPUTS:
%   avgEEG - average EEG values  (nSamples x 1 vector)
%   errEEG - error of average EEG values (std error of mean)
%   EEGts - timetamps of each 'bin' (sample point) in the PETH
%   EEGsamples - individual EEG snippets around each event before averaging (length(T) x nSamples array) 
%
% Peri-Event-Time Histogram (PETH) of an EEG file with respect to a list of trigger-event timestamps T.
% The window size around T (at 0) can be asymmetric and ranges from [T - tlow_msec, T + thigh_msec].
% The number of bins depends on the sample frequency in the EEG file and the length of your window.
%
% PL 2001, last modified '02 by MN


eeg = {};
eeg_ts = {};

nT = length(T);
ix0 = zeros(nT,1);
maxinfimum = Inf;
minsupremum = Inf;

for i = 1:nT
    [eeg_ts{i}, eeg{i}, sFreq] = ReadCR_nt_partial_load( eeg_fname, T(i)-abs(tlow_msec)*10, T(i)+abs(thigh_msec)*10);
    
    % interplate timestamps of each block
    blockSize = 512;
    nBlocks = size(eeg{i},1);
    eeg{i} = reshape(eeg{i}',blockSize*nBlocks,1);
    dT = 10000/sFreq;    % in tstamps
    
    ts = eeg_ts{i};
    ts = [ts;ts(end) + 512*dT];     
    eeg_ts{i} = zeros(size(eeg{i}));
    for iBlock = 1:(nBlocks)
        eeg_ts{i}((blockSize * (iBlock-1) + 1):(blockSize * iBlock)) = ...
            linspace(ts(iBlock), ts(iBlock+1) - dT, blockSize);
    end
    
    [dtsmin, ix0(i)] = min(abs(eeg_ts{i}-T(i) ));    % find index ix0 of closest timestamp to 0 
    if ix0(i)-1 < maxinfimum
        maxinfimum = ix0(i)-1;
    end
    if (length(eeg_ts{i}) - ix0(i)) < minsupremum
        minsupremum = (length(eeg_ts{i}) - ix0(i));
    end
    
end

% align all nT EEG chunks (of varying lengths) to ix0
nSamples = maxinfimum + minsupremum + 1;
EEGts = eeg_ts{1}(ix0(1)-maxinfimum : ix0(1)+minsupremum)-T(1);
if(length(EEGts) ~= nSamples)
    error([ num2str(length(EEGts)) ' ~= ' num2str(nSamples) '!!']);
end

sum_EEGsamples = zeros(nSamples,1);
EEGsamples={};
for i = 1:nT
    EEGsamples{i} = eeg{i}(ix0(i)-maxinfimum : ix0(i)+minsupremum);
    sum_EEGsamples = sum_EEGsamples + EEGsamples{i};
end

avgEEG = sum_EEGsamples/nT;

SS_EEGsamples = zeros(nSamples,1);
for i = 1:nT
    SS_EEGsamples = SS_EEGsamples + (EEGsamples{i}-avgEEG).^2;
end

if nT > 1
    errEEG = sqrt(SS_EEGsamples/(nT-1))/sqrt(nT);
else
    errEEG = [];
    disp('EEG standard errors not defined- only one event to average over!');
end