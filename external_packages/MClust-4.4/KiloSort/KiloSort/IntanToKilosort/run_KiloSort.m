function run_KiloSort(fpathIn, fpath, data_fname, chanMap_fname, config_fname)
% Run Kilosort on raw intan data.  Works for both tetrodes and Silicon
% probes as long as you use the correct channel map file.
%
% USAGE
%   run_KiloSort(fpathIn, fpath, data_fname, chanMap_fname, config_fname);
%
% INPUT
%   fpathIn - directory path containing raw Intan files
%   fpath - directory to generate output kilosort data into
%   data_fname - filename for the single Kilosort-format binary data
%   chanMap_fname - filename of channel map
%   config_fname - filename for the Kilosort configuration file
%
% OUTPUT
%   Runs kilosort which generates the rez.mat file and python files which
%   can be read by phy.
%
% EXAMPLE
%   run_KiloSort('C:\Data\R422\R422-2017-05-18', ...
%       'R422-2017-05-18-KilosortRaw-TTs', ...
%       'C:\...\KiloSort\channel_map_24TT.txt', ...
%       'C:\...\KiloSort\config_KiloSort_24TT.m');
%
% June 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities

% Paths and Filenames
%fpath = 'C:\Users\hasz\Documents\Data-temp\Kilosort-32Si-test\'; % where on disk do you want the simulation? ideally and SSD...
%data_fname = 'R422-2017-05-18-32SiKilosortRaw.dat'; %path to Kilosort-formatted data
%chanMap_fname = 'C:\Users\hasz\Documents\MATLAB\BMHcodeset\Neural\KiloSort\channel_map_32Si.txt'; %filename of channel map
%config_fname = 'C:\Users\hasz\Documents\MATLAB\BMHcodeset\Neural\KiloSort\config_KiloSort.m'; % configuration filename

% Check if there is existing Kilosort data in the output folder
KO = cell(5,1); %cell array of Kilosort Output files
KO{1} = fullfile(fpath, 'spike_times.npy');
KO{2} = fullfile(fpath, 'spike_clusters.npy');
KO{3} = fullfile(fpath, 'cluster_groups.csv');
KO{4} = fullfile(fpath, 'spike_templates.npy');
KO{5} = fullfile(fpath, 'templates.npy');
for iF = 1:length(KO)
    if exist(KO{iF}, 'file')
        warning([KO{iF} ' already exists. Deleting it.'])
        delete(KO{iF});
    end
end

% Check if data_fname is full path
if ~contains(data_fname, '\') %not full path
    data_fname = fullfile(fpath, data_fname);
end

% Convert raw data to Kilosort format
ConvertRawToKiloSort(fpathIn, data_fname, chanMap_fname)

% Run the configuration file, it builds the structure of options (ops)
run(config_fname)
ops.fbinary = data_fname; %set the data filename

% This part makes the channel map for Kilosort
chanMapOut = fullfile(fpath, 'chanMap.mat');
WriteKiloSortChannelMap(chanMap_fname, chanMapOut);

% Run KiloSort
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% Save python results file for Phy
rezToPhy(rez, fpath);
save(fullfile(fpath,  'rez.mat'), 'rez');

% now fire up Phy and check these results.
% Open an Anaconda prompt, navigate to fpath, and run:
%   activate phy
%   phy template-gui params.py

end