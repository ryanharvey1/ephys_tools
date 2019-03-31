function ConvertRawToKiloSort(dirIn, fnameOut, chanMapFname)
% Convert raw Intan .dat files from each channel into a single 
% KiloSort .dat file.
%
% USAGE
%   ConvertRawToKiloSort(dirIn, fnameOut, chanMapFname);
%
% INPUT
%   dirIn - directory path containing raw Intan files
%   fnameOut - filename for single output Kilosort-format binary data
%   chanMapFname - filename of channel map
%
% OUTPUT
%   Writes raw data in KiloSort-readable format to fnameOut.
%
% June 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities

% Settings
ChunkSize = 1000000; %number of samples per chunk

% Load Intan -> TT/Si Probe channel map
eval(fileread(chanMapFname));

% Get number of samples in entire recording
fileName = [dirIn filesep 'amp-' port(1) ...
    sprintf('-%0.3d.dat', channel_map(1))];
sz = dir(fileName);
Nts = sz.bytes/2; %total number of int16 samples in each raw file

% Open raw Intan files
Nc = numel(channel_map); %number of channels
fid_in = cell(Nc, 1);  %cell array of file pointers
for iC = 1:Nc %for each channel,
    fname_in = [dirIn filesep 'amp-' port(iC) ...
        sprintf('-%0.3d.dat', channel_map(iC))];
    fid_in{iC} = fopen(fname_in, 'r'); %open the raw file
end

% Open Kilosort file to write
fid_out = fopen(fnameOut, 'a'); 

% Write Intan .dat files to single KiloSort .dat file in chunks
S = zeros(Nc, ChunkSize, 'int16'); %array to store chunks
SamplesRead = 0;
iChunk = 1;
Nchunk = ceil(Nts/ChunkSize);
fprintf('Saving Raw Intan data to Kilosort Format ... 000/%0.3d', Nchunk);
while SamplesRead<Nts %while we haven't processed all samples in the files,

    fprintf('\b\b\b\b\b\b\b%0.3d/%0.3d', iChunk, Nchunk);
    
    % Read in this chunk
    Ntc = min(ChunkSize, Nts-SamplesRead); % Number samples for this chunk 
    for iC = 1:Nc % read a chunk from each channel
        S(iC,1:Ntc) = fread(fid_in{iC}, Ntc, 'int16');
    end
    
    % Write to the single KiloSort-formatted binary file
    fwrite(fid_out, S(:,1:Ntc), 'int16');
    
    % Update how many samples we've read
    SamplesRead = SamplesRead + Ntc;
    iChunk = iChunk + 1;
    
end

% Close files
for iC = 1:Nc
    fclose(fid_in{iC});
end
fclose(fid_out);
fprintf('\nDone.\n');

end