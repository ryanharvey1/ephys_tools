function [Ns, Nclu] = WriteSingleTetSpikesFromKilosort(ST, oC, rawFname, AT, fpath, chanMap_fname, ssn, TTN)
% Write a .spikes file for a single tetrode from kilosort data, and a
% corresponding .clu file for that tetrode based on kilosort clusters.
% Before using this function you have to compile SpikeFilter! (see below)
%
% USAGE
%   WriteSingleTetSpikesFromKilosort(ST, C, ...
%       fpathIn, fpath, chanMap_fname, ssn)
%
% INPUT
%   ST - spike times (in sample #s) (Nspikes-length vector)
%   oC - cluster IDs for each spike (Nspikes-length vector)
%   rawFname - filename of raw data file which was fed to Kilosort
%   AT - aligned time (time of each sample, aligned from timestamps)
%   fpath - directory path for output .spikes and .clu files
%   chanMap_fname - filename of channel map
%   ssn - R###-YYYY-MM-DD string
%   TTN - tetrode number for this tetrode
%
% OUTPUT
%   Writes a .spikes file with the waveforms for each spike on a tetrode,
%   and a .clu file for each tetrode with the cluster IDs.
%
% NOTE
%   Before using this function, you must compile SpikeFilter.c with:
%   
%   mex SpikeFilter.c
%
% ANOTHER NOTE:
%   This takes kind of a long time, and I'm not sure why...  ~400s for a
%   single TT.  It may be the common median referencing.
%
% July 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities

%% Settings
CMR = true;   % whether to perform common median referencing
locut = 600;  % lower freq for bandpass
hicut = 7500; % upper freq for bandpass
window = 64;  % window size
offset = 20;  % offset from beginning of window to spike time
ChunkSize = 1000000; % number of samples per chunk


%% Load and process raw data

% Load Intan -> TT/Si Probe channel map
eval(fileread(chanMap_fname));

% Get number of channels and samples in entire recording
Nc = numel(channel_map); %number of channels
sz = dir(rawFname);
Nts = sz.bytes/2/Nc; %total number of int16 samples in each raw file

% Open Kilosort file to write
fid = fopen(rawFname, 'r'); 

% Load data from desired TTs from raw KiloSort .dat file
% and perform common median referencing
Ncps = size(channel_map,2); %number of channels per sensor/TT
S = zeros(Ncps, Nts, 'int16'); %allocate array for this TT
SamplesRead = 0;
while SamplesRead<Nts %while we haven't processed all samples in the file,
    
    % Read in this chunk
    Ntc = min(ChunkSize, Nts-SamplesRead); % Number samples for this chunk 
    Sr = fread(fid, [Nc Ntc], 'int16=>int16'); %read this chunk from file
    
    % Copy data for desired TT
    S(:, SamplesRead+1:SamplesRead+Ntc) = ...
        Sr(find(shank==TTN), :); %#ok<FNDSB>
    
    % Subtract median
    if CMR
        S(:, SamplesRead+1:SamplesRead+Ntc) = ...
            S(:, SamplesRead+1:SamplesRead+Ntc) - median(Sr);
    end
        
    % Update how many samples we've read
    SamplesRead = SamplesRead + Ntc;
    
end

% Close raw file
fclose(fid);

% Transpose S (such that each column is data from 1 TT)
S = S';

% Bandpass
Sf = int16(SpikeFilter(double(S), locut, hicut, fs));
clear S


%% Write .spikes file

% Don't use spikes which are w/i window of beginning or end
keep = ST>window & ST<(Nts-window);
ST = ST(keep);
oC = oC(keep);

%Extract spikes @ times from Kilosort
Ns = length(ST); %number of spikes
spikes = zeros(window, Ncps, Ns, 'int16');
for iS = 1:Ns %for each spike,
    spikes(:,:,iS) = Sf(ST(iS)-offset+1:ST(iS)-offset+window, :);
end

% Write the .spikes file
spikesFname = fullfile(fpath, [ssn sprintf('-TT%0.2d.spikes',TTN)]);
fid = fopen(spikesFname, 'w'); 
fwrite(fid, uint32(Ns), 'uint32');      % Number of spikes
fwrite(fid, uint16(Ncps), 'uint16');    % Number of channels
fwrite(fid, uint16(window), 'uint16');  % Length of spike waveform window
fwrite(fid, AT(ST+1), 'double');        % Spike times
fwrite(fid, spikes, 'int16');           % Spike waveforms
fclose(fid);


%% Write .clu file

% Get unique cluster IDs for this TT
C = nan(size(oC));
uC = unique(oC); %unique cluster IDs
Nclu = length(uC); %number of clusters
for iC = 1:Nclu %for each unique cluster
    C(oC==uC(iC)) = iC; %give new monotonically increasing ID
end
C = [Nclu; C(:)]; %first line is number of classes found

% Write the .clu file
cluFname = fullfile(fpath, [ssn sprintf('-TT%0.2d.clu',TTN)]);
dlmwrite(cluFname, C, '%d');


end