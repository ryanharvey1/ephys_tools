function SaveSpikeTimesFromKilosort(timePath, fpathIn, fpath, sensor)
% Write .t files (containing spike times) for each kilosort spike cluster
% Before running this function, you should run kilosort on your raw data 
% (run_KiloSort.m), and sort/merge/split the clusters in Phy.
%
% USAGE
%   SaveSpikeTimesFromKilosort(timePath, fpath, chanMapFname)
%
% INPUT
%   timePath - path containing timestamps (board-DIN-00.dat etc)
%   fpathIn - directory containing KiloSort output .npy files
%   fpath - directory path for output .t files
%   sensor - string with the sensor (e.g. "24TT" or "SiPL" or "SiIL")
%
% OUTPUT
%   Writes a .t file (with spike times) for each spike cluster.
%
% EXAMPLE
%   SaveSpikeTimesFromKilosort(...
%       'F:\Data\Raw\R422\R422-2017-05-18', ...
%       'F:\Data\Processed\R422\R422-2017-05-18', ...
%       'C:\Users\hasz\Documents\MATLAB\KiloSort\channel_map_24TT.txt');
%
% August 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities


%% Filenames

% .npy files
spike_time_fname = fullfile(fpathIn, 'spike_times.npy');
spike_clusters_fname = fullfile(fpathIn, 'spike_clusters.npy');
cluster_groups_fname = fullfile(fpathIn, 'cluster_groups.csv');

% SSN
ind = strfind(fpath,filesep);
ind = ind(end)+1;
ssn = fpath(ind:ind+14);


%% Load spike times
ST = readNPY(spike_time_fname);


%% Load clusters

% Read cluster IDs for each spike
C = readNPY(spike_clusters_fname);

% Read cluster info
cid = [];
cinf = [];
fid = fopen(cluster_groups_fname);
fgetl(fid); %discard 1st line
tline = fgetl(fid);
while ischar(tline) %until we've reached the end of file,
    twords = textscan(tline, '%s');
    cid = [cid; str2double(twords{1}{1})]; %#ok<AGROW>
    if strcmp('noise', twords{1}{2}) %noise
        cinf = [cinf; 0];  %#ok<AGROW>
    elseif strcmp('mua', twords{1}{2}) %mua
        cinf = [cinf; 1];  %#ok<AGROW>
    elseif strcmp('good', twords{1}{2}) %good
        cinf = [cinf; 2];  %#ok<AGROW>
    else
        cinf = [cinf; 3];  %#ok<AGROW>
    end
    tline = fgetl(fid); %get next line
end
fclose(fid); %close file


%% Synchronize Intan and Matlab time

% Get aligned time vector from timestamps
clkFname = [timePath filesep 'board-DIN-00.dat'];
dataFname = [timePath filesep 'board-DIN-01.dat'];
if exist(clkFname, 'file')==2 && exist(dataFname, 'file')==2 %timestamps?
    T = GenerateSyncedNeuralTimes(clkFname, dataFname); % get time from em
else %otherwise just use default Intan time
    Tname = [timePath filesep 'time.dat']; %filename for time
    fid = fopen(Tname, 'rb');
    sz = dir(Tname);
    T = fread(fid, sz.bytes/4, 'int32');
    fclose(fid); %close the file
    T = double(T);
end


%% Write a .t file for each good cluster

iG = 1; %number of good units
Ns = 0;
for iC = 1:length(cid) %for each cluster,
    if cinf(iC)==2 % only keep good clusters
        tST = T(ST(C==cid(iC))); % times of spikes in this cluster
        fname = sprintf('%s\\%s-%s_%0.2d.t', fpath, ssn, sensor, iG); % .t filename
        WriteSpikeTimesFile(fname, tST); % write the .t file
        iG = iG + 1; %add to good clusters
        Ns = Ns + length(tST); %count total number of spikes
    end
end
fprintf('  Wrote .t files for %d cells ( %d spikes )\n', iG-1, Ns);


end
