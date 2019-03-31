function ExtractSpikesFromKilosort(rawFname, timePath, fpathIn, fpath, chanMapFname)
% Write .spikes files (containing waveforms and spike times) and .clu files
% (containing cluster IDs for each spike) for all tetrodes from kilosort 
% data and raw Intan .dat files.  Before running this function, you should
% run kilosort on your raw data (run_KiloSort.m), and sort the clusters
% into good or noise in Phy.
%
% USAGE
%   ExtractSpikesFromKilosort(fpathIn, fpath, chanMapFname)
%
% INPUT
%   rawFname - filename of raw data file which was fed to Kilosort
%   timePath - path containing timestamps (board-DIN-00.dat etc)
%   fpathIn - directory containing the KiloSort output .npy files
%   fpath - directory path for output .spikes and .clu files
%   chanMapFname - filename of channel map
%
% OUTPUT
%   Writes .spikes and .clu files for each tetrode.
%
% EXAMPLE
%   ExtractSpikesFromKilosort(...
%       'F:\Data\Raw\R422\R422-2017-05-18\R422-2017-05-18-KilosortRaw-24TT.dat', ...
%       'F:\Data\Raw\R422\R422-2017-05-18', ...
%       'F:\Data\Processed\R422\R422-2017-05-18', ...
%       'C:\Users\hasz\Documents\MATLAB\KiloSort\channel_map_24TT.txt');
%
% July 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities


%% Filenames

% .npy files
spike_time_fname = fullfile(fpathIn, 'spike_times.npy');
spike_clusters_fname = fullfile(fpathIn, 'spike_clusters.npy');
cluster_groups_fname = fullfile(fpathIn, 'cluster_groups.csv');
spike_templates_fname = fullfile(fpathIn, 'spike_templates.npy');
spike_template_dat_fname = fullfile(fpathIn, 'templates.npy');

% SSN
ind = strfind(fpath,filesep);
ind = ind(end)+1;
ssn = fpath(ind:ind+14);

%% Load Intan -> TT/Si Probe channel map
eval(fileread(chanMapFname));

%% Load spike times
ST = readNPY(spike_time_fname);

%% Find what TT each spike is on from templates

% Load templates for each spike
Temps = readNPY(spike_templates_fname);

% Find what TT each template is on from spike templates
TTids = shank(connected); % TT # of each channel
TP = readNPY(spike_template_dat_fname); %spike templates (Ntemplates-by-Nsamples-by-Nchannels matrix)
tempTT = nan(size(TP,1),1); %what TT each template is on
for iT = 1:size(TP,1) %for each spike template,
    enC = squeeze(sum(power(TP(iT,:,:),2),2)); %energy of this template on each channel
    [~,mind] = max(enC); %channel w/ max energy
    tempTT(iT) = TTids(mind); %tetrode ID for this template
end

% Find what TT each spike was on
TTs = Temps;
uTemps = unique(Temps); %unique template IDs
for iT = 1:length(uTemps)
    TTs(Temps==uTemps(iT)) = tempTT(uTemps(iT)+1);
end

%% Find spikes to keep based on cluster IDs

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

% Find out what spikes to keep
keep = nan(size(C));
for iC = 1:length(cid) %for each cluster,
    if cinf(iC)==2 %only keep good clusters
        keep(C==cid(iC)) = true;
    else
        keep(C==cid(iC)) = false;
    end
end


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


%% Write a .spikes and .clu file for each TT

tets = unique(shank);
t1 = tic;
totu = 0;
for iT = 1:length(tets)
    fprintf('TT %d / %d ...\n', iT, length(tets));
    t2 = tic;
    tinds = keep & TTs==tets(iT);
    if nnz(tinds)>0
        tT = ST(tinds); %spike times for spikes on this TT
        tC = C(tinds);  %cluster IDs for spikes on this TT
        [Ns, Nclu] = WriteSingleTetSpikesFromKilosort(tT, tC, ... %write 
            rawFname, T, fpath, chanMapFname, ssn, tets(iT));     %spikes
        fprintf('  Extracted %d spikes\n  %d putative units\n', Ns, Nclu);
        disp(['  Elapsed time: ' num2str(toc(t2)) ' seconds.']);
        totu = totu + Nclu; %count total number of units
    else
        fprintf('  No spikes detected on TT %d\n', iT);
    end
end
disp(['Total elapsed time: ' num2str(toc(t1)) ' seconds.']);
disp(['Total putative units: ' num2str(totu)]);

end
