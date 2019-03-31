%testingPhasePrecession 

% Run postprocesMClust_v9.m to   
%
%         % Compute Phase Lock
%         Stats=EEGWorking2(eegfile{1},spks_VEL,StartofRec,EndofRec,event);
%
% Then play through sections

% Load EEG data:
addpath(genpath('F:\Users\BClarkLab\GoogleDrive\MatlabDir\NSMA_Toolbox'));
eeg_tsd=ReadCR_tsd(eegfile{1});
%% Filter for theta frequency:
eeg_theta = Filter4Theta(eeg_tsd,6,11);
% Find theta phase during each spike:
[ph, thpeaks] = ThetaPhase(S, eeg_theta,[],[]);
% Plot theta phase versus rat position for all spikes (of cell# c)
c=1;
phase=Data(ph{c});
plot(Vhs_spikes{c}, phase(:,1), 'k.')
hold on; plot(Vhs_spikes{c}, phase(:,1)+1, 'b.')



%%
figure(5)
lonWrapped = wrapTo360(rad2deg(spks_VEL_working(:,2)));
plot(spks_VEL(:,2), lonWrapped, 'k.')
hold on; plot(spks_VEL(:,2), lonWrapped+360, 'r.')
ylim([0 720])

figure(7)
% lonWrapped = wrapTo360(rad2deg(spks_VEL_working(:,2)));
plot(spks_VEL(:,2), spks_VEL_working(:,2), 'k.')
hold on; plot(spks_VEL(:,2), spks_VEL_working(:,2)+1, 'r.')

