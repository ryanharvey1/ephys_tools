
addpath(genpath('Z:\Matlab2\User_Specific\cowen'))
addpath(genpath('u:\matlab'))
addpath('U:\matlab')

data_dir = pwd;
ts_txtfile = [data_dir '\epochs_ts.txt'];

f = findstr(data_dir, '200');
g = findstr(data_dir, 'Dorsal_ventral\');
h= findstr(data_dir, 'Datasets\');
dmv_str = data_dir(h+9);
rat = data_dir(g+15:g+18);
sess_num = data_dir(f:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading spikes- Need to find the tt with the most cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InterneuronPyramidalCellIndex = []; %Documents which cells are which: PyramidalCells == 1; Interneurons == 0
data_dir = pwd;
tfile_dir = [data_dir '\Pyramidal\Good'];
get_spikes;

tetrodenames = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
CompleteBestTT = [];
Nfiles = [];
for iTT = 1:length(tetrodenames);
    Nfiles = [];
    for ii = 1:length(tfiles)
        foundtheTT = findstr(tetrodenames{iTT}, tfiles{ii});
        if ~isempty(foundtheTT)
            Nfiles = [Nfiles; iTT];
        end
    end
    if ~isempty(Nfiles);
        if ~isempty(CompleteBestTT)
            if length(Nfiles) > CompleteBestTT(:,2)
                CompleteBestTT = [unique(Nfiles), length(Nfiles)];
            end
        else
            CompleteBestTT = [unique(Nfiles), length(Nfiles)];
        end
    end
end
cd(data_dir)
CSCtoFind = ['CSC0' num2str(CompleteBestTT(:,1)) '.ncs'];
if length(CSCtoFind) > 9;
    CSCtoFind = ['CSC' num2str(CompleteBestTT(:,1)) '.ncs'];
end
eegfiles=findfiles(CSCtoFind);
if isempty(eegfiles)
    CSCtoFind = ['CSC' num2str(CompleteBestTT(:,1)) '.ncs'];
    eegfiles=findfiles(CSCtoFind);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ref_theta_csc = 'CSC8.ncs';
ref_theta_csc = CSCtoFind;

pushdir(data_dir);
%----------------------------------------
% load epochs
[t, mz_config] = load_epoch_ts(ts_txtfile);

switch mz_config
    case 0      % small track, big track
        t_mz{1} = t(7:8);
        t_mz{2} = t(3:4);
    case 1      % big track, small track
        t_mz{1} = t(3:4);
        t_mz{2} = t(7:8);
end


track_str = {'BigTrack', 'SmallTrack'};
%-------------------------------------------------------
% process EEG file
iTrack = 1;   % loop over big=1 and small=2 tracks

tic;
% load smoothed (1sec smoothing) angular tracker data from .mat file
switch iTrack
    case 1
        load([data_dir '\spike_rast_vars_big.mat']);        % loads Xr_ts and thPos{1:2} angular position
        posTh = thPos{1};
        track_diam = 121.5;
    case 2
        load([data_dir '\spike_rast_vars_small.mat']);
        posTh = thPos{2};
        track_diam = 53.3;
end
clear thPos;

% estimate instantanous velocity from smoothed angular position
dts = diff(Xr_ts);
vel_ts = Xr_ts(1:end-1)+dts/2;      % set velocity timestamps at midpoints between two video frames
vel_cm = track_diam/2 * abs(diff(posTh))./(dts/10000);  %  instantanous |velocity| in cm/sec
vel_cm(vel_cm > 100) = 100;     % truncate extreme velocity outliers to 100 cm/sec

% get eeg data from theta reference csc.dat files and
% filter each ref eeg for theta
disp(['Reading CSC file ... ' ref_theta_csc]);
eeg_header = Read_nlx_header(ref_theta_csc);
iTrack = 2;
lengthofepoch = (t_mz{iTrack}(2).*100) - (t_mz{iTrack}(1).*100);
increment = lengthofepoch./100;
%timestorunrand = (t_mz{iTrack}(1).*100):10000:(t_mz{iTrack}(2).*100);
mega_mtx = [];

for i = 0:20;
    i
    if ((t_mz{iTrack}(1).*100)+(i.*increment) + increment) <=  (t_mz{iTrack}(2).*100);
        clear T EEG
        [T, EEG] = Read_CR_files(ref_theta_csc, 2000, [(t_mz{iTrack}(1).*100)+ increment.*i,(t_mz{iTrack}(1).*100)+ increment.*(i+1)]);
        [TFR,timeVec,freqVec] = traces(EEG,[1:1:200],2000,9);
        mega_mtx = [mega_mtx; squeeze(mean(TFR'))];
    end
end

figure;
semilogy(1:200, mean(mega_mtx))
hold on;            semilogy(1:200,mean(mega_mtx) + Sem(mega_mtx),'r')
hold on;            semilogy(1:200,mean(mega_mtx) - Sem(mega_mtx),'r')

figure;
hold on;
semilogy(1:200, mean(mega_mtx)./sum(mean(mega_mtx)), 'r')


% figure; errorbar([1:1:200], mean(booted_mtx), std(booted_mtx)./sqrt(size(booted_mtx,1)))
% norm_M = Z_scores(booted_mtx);
% 
% new_mtx = [];
% for ii = 1:size(booted_mtx,2);
%     new_mtx = [new_mtx, booted_mtx(:,ii)./sum(booted_mtx(:,ii))];
% end
% figure; errorbar([1:1:200], mean(new_mtx), std(new_mtx)./sqrt(size(new_mtx,1)))
% 
% 
% figure; bar(sum(TFR'))
% figure; plot(log10(squeeze(mean(TFR'))))
% %for ii  =1:size(TFR,2)