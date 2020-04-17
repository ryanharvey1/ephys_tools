function ops = StandardConfig_KS2Wrapper(basepath)
%% Standard configuration for kilosort2 wrapper
%  Remember that this configuration was optimized for striatal recordings,
%  so make sure to adapt to your recordings
%
% Eliezyer de Oliveira 2020/03/02


%directory to write temporary file (should be a SSD for speed)
rootH = basepath; 

%loading xml file 
d   = dir('*.xml');
if ~isempty(d)
    par = LoadXml(fullfile(basepath,d(1).name));
else
    error('the .xml file is missing')
end

%range to spike sort
ops.trange = [0 Inf];

ops.NchanTOT  = par.nChannels; % total number of channels in your recording
%channel map file location
ops.chanMap   = fullfile(basepath,'chanMap.mat');
%sampling rate of the .dat file
ops.fs        = par.SampleRate;  
% frequency for high pass filtering (Hz)
ops.fshigh    = 500;
% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0; 
% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [10 4];  
% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
ops.lam = 5;  
% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9; 
% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = -Inf; %-Inf works? striatum has pretty low firing rates, that's why I'm using this
% number of samples to average over (annealed from first to second value) 
ops.momentum = [20 400]; 

% spatial constant in um for computing residual variance of spike
ops.sigmaMask = 30; 

ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; 
%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction. 
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
ops.whiteningRange      = 32; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

%%