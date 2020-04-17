 function KS2Wrapper(basepath)
%Kilosort 2 wrapper for neuroscope users. It extract informations from your
%.xml to perform spike sorting.
%
%Eliezyer de Oliveira 02/03/2020


%% some flags for the script
neurosuite_files = 0;

%% you need to change most of the paths in this block
if nargin<1
    basepath = pwd;
end

%loading configuration structure
ops = StandardConfig_KS2Wrapper(basepath);

%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', basepath)

% is there a channel map file in this folder?
fs = dir(fullfile(basepath, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(basepath, fs(1).name);
else
    createChannelMapFile_Local(basepath)
end
%IF NOT, CREATE A CHANNEL MAP BASED ON THE XML

% find the binary file
fs          = [dir(fullfile(basepath, '*.bin')) dir(fullfile(basepath, '*.dat'))];
ops.fbinary = fullfile(basepath, fs(1).name);

% preprocess data to create temp_wh.dat
disp('Preprocessing data')
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
disp('drift correction...')
rez = clusterSingleBatches(rez);

%creating a folder to save kilosort results every time
temp = ['Kilosort2_' datestr(clock,'yyyy-mm-dd_HHMMSS')];
ks2_folder = fullfile(basepath, temp);
mkdir(ks2_folder);

% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(ks2_folder, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1); %remove later?

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

%% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, ks2_folder);

%correcting params.py file to access the .dat file
fid = fopen(fullfile(ks2_folder,'params.py'), 'w');
[~, fname, ext] = fileparts(rez.ops.fbinary);
fprintf(fid,['dat_path = ''../',fname ext '''\n']);
fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
fprintf(fid,'dtype = ''int16''\n');
fprintf(fid,'offset = 0\n');
if mod(rez.ops.fs,1)
    fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
else
    fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
end
fprintf(fid,'hp_filtered = False');
fclose(fid);

%% saving results to matlab

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(ks2_folder, 'rez2.mat');
save(fname, 'rez', '-v7.3');

%% creating neurosuite files if specified in the variable neurosuite_files
 if neurosuite_files
    ConvertKilosort2Neurosuite_KSW(rez) 
 end