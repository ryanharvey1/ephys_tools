
function run_ks2(basepath)

%% you need to change most of the paths in this block
rootZ = basepath; % the raw data binary file is in this folder
rootH = basepath; % path to temporary binary file (same size as data, should be on fast SSD)
pathToYourConfigFile = 'D:\Users\BClarkLab\ephys_tools\external_packages\Kilosort2'; % take from Github folder and put it somewhere else (together with the master_file)

ops.trange    = [0 Inf]; % time range to sort
ops.NchanTOT  = 64; % total number of channels in your recording

run(fullfile(pathToYourConfigFile, 'StandardConfig_E14x16.m'))
ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD


%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% find the binary file
fs          = dir(fullfile(rootZ, '*.dat'));
ops.fbinary = fullfile(rootZ, fs(1).name);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);
%

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);
                 
% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% remove double-counted spikes - solves issue in which individual spikes are assigned to multiple templates.
% See issue 29: https://github.com/MouseLand/Kilosort/issues/29
rez = remove_ks2_duplicate_spikes(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% decide on cutoff
rez = set_cutoff(rez);

% eliminate widely spread waveforms (likely noise)
rez.good = get_good_units(rez);

fprintf('found %d good units \n', sum(rez.good>0))

%creating a folder to save kilosort results every time
temp = ['Kilosort2_' datestr(clock,'yyyy-mm-dd_HHMMSS')];
ks2_folder = fullfile(basepath, temp);
mkdir(ks2_folder);

% write to Phy
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

%% if you want to save the results to a Matlab file...

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% final time sorting of spikes, for apps that use st3 directly
[~, isort]   = sortrows(rez.st3);
rez.st3      = rez.st3(isort, :);

% Ensure all GPU arrays are transferred to CPU side before saving to .mat
rez_fields = fieldnames(rez);
for i = 1:numel(rez_fields)
    field_name = rez_fields{i};
    if(isa(rez.(field_name), 'gpuArray'))
        rez.(field_name) = gather(rez.(field_name));
    end
end

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(ks2_folder, 'rez2.mat');
save(fname, 'rez', '-v7.3');


end