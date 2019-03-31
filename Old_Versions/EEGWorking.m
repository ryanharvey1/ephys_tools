function [Stats]=EEGWorking(eegfile,spks_VEL,StartofRec,EndofRec,event)
% Imports eegfile dirs & timestamps containing spike values
% Exports Mean Resultant Vector, Rayleigh's Test (p,zval) for delta-high gamma

% Ryan Harvey FEB 2017

addpath(genpath( 'd:\Users\BClarkLab\GoogleDrive\MatlabDir\eeglab12_0_1_0b'));
addpath(genpath( 'F:\Users\BClarkLab\GoogleDrive\MatlabDir\eeglab12_0_1_0b'));

% extract CSC data with mex
if ismac==1
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples]= Nlx2MatCSC_v3(eegfile, [1 1 1 1 1], 0, 1);
else
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples]= Nlx2MatCSC(eegfile, [1 1 1 1 1], 0, 1);
end
% reshape Samples variable into a vector
EEGreshaped = reshape(Samples,1,length(Samples(:)));

% increase the number of timestamps by a factor of 512 so that there is a timestamps for every data point
factor = 512;
EEGtsUpSampled = zeros(factor,length(Timestamps));
EEGtsUpSampled(1,:)=Timestamps;
TSdiff=diff(Timestamps);
TSdiff(end+1)=TSdiff(end);
TSdiffIncrements=TSdiff/factor;
for i=1:length(Timestamps)
    for j=2:factor
        EEGtsUpSampled(j,i)=Timestamps(i)+TSdiffIncrements(i)*(j-1);
    end
end
EEGtsUpSampled = (EEGtsUpSampled(:))';

% downsample to 1000 Hz (is a little nicer to work with than 1523 Hz (Taube lab) or 2034 Hz (McN lab))
NewSFreq = 1000;

% calculate number of samples in down sampled data
DurationData = (EEGtsUpSampled(end)-EEGtsUpSampled(1))/1000000; % in seconds
NumNSamples = floor(DurationData*NewSFreq);

% create new timestamps
index = 1:NumNSamples;
EEG_DownSampledTimestamps(index) = EEGtsUpSampled(1)+(1000000*1/NewSFreq)*(index-1);

% interpolate data points
EEG_DownSampledData = interp1(EEGtsUpSampled, EEGreshaped, EEG_DownSampledTimestamps);

% RESTRICT DATA BY START AND END OF EVENT
if (length(StartofRec) + length(EndofRec))>2
    EEG_DownSampledData=EEG_DownSampledData(EEG_DownSampledTimestamps>StartofRec(event) & EEG_DownSampledTimestamps<EndofRec(event));
    EEG_DownSampledTimestamps=EEG_DownSampledTimestamps(EEG_DownSampledTimestamps>StartofRec(event) & EEG_DownSampledTimestamps<EndofRec(event));
end

% PHASE LOCKING
for freq=1:6
    % filter between 4 and 12 Hz.
    if freq==1; EEGthetaData = eegfilt(EEG_DownSampledData,NewSFreq,1,4);  % DELTA
    elseif freq==2; EEGthetaData = eegfilt(EEG_DownSampledData,NewSFreq,4,8);  % THETA
    elseif freq==3; EEGthetaData = eegfilt(EEG_DownSampledData,NewSFreq,8,13);  % ALPHA
    elseif freq==4; EEGthetaData = eegfilt(EEG_DownSampledData,NewSFreq,13,30);  % BETA
    elseif freq==5; EEGthetaData = eegfilt(EEG_DownSampledData,NewSFreq,30,150);  % GAMMA
    elseif freq==6; EEGthetaData = eegfilt(EEG_DownSampledData,NewSFreq,150,200); % HIGH GAMMA
    end
    % determine theta phase using the Hilbert transform
    ThetaPh = angle(hilbert(-EEGthetaData));
    
    ts_spike_working = interp1(EEG_DownSampledTimestamps', EEG_DownSampledTimestamps', spks_VEL(:,1), 'nearest');
    % ASSIGN SPIKES TO TIMESTAMPS
    [~,ia,~] = intersect(EEG_DownSampledTimestamps', ts_spike_working);
    N = size(EEG_DownSampledTimestamps'); spikeVec = zeros(N(1,1),1); spikeVec(ia) = 1;
    % CONCAT TIMESTAMPS, PHASE, AND SPIKE [0 or 1]
    EEG_DownSampledTimestamps_working = [EEG_DownSampledTimestamps' ThetaPh' spikeVec(:,1)];
    % ISOLATE SPIKING
    spks_VEL_working = EEG_DownSampledTimestamps_working(EEG_DownSampledTimestamps_working(:,3) == 1,:);
        
    
    % CALCULATE SPIKES OF PHASE STATS
    % R LENGTH
    Stats.MeanResultantPhase(freq,:)=circ_r(spks_VEL_working(:,2)); 
    % RAYLEIGH TEST
    [Stats.RayleighsTest.PVal(freq),Stats.RayleighsTest.ZValue(freq)] = circ_rtest(spks_VEL_working(:,2)); 
    % CONFIDENCE INTERVAL
    [~,Stats.Boundaries(freq,:)] = CircularConfidenceIntervals(spks_VEL_working(:,2)); 
 
    % VARIANCE
    [Variance,StandardDev] = CircularVariance(spks_VEL_working(:,2)); 
    Stats.Variance(freq)=rad2deg(Variance)*2;
    Stats.StandardDev(freq)=rad2deg(StandardDev);
    % MEAN
    stats = circ_stats(spks_VEL_working(:,2));
    if rad2deg(stats.mean)<0;
        Stats.Meanz(freq)=rad2deg(stats.mean)+360;
    elseif rad2deg(stats.mean)>0;
        Stats.Meanz(freq)=rad2deg(stats.mean);
    end
    % MEDIAN
    if rad2deg(stats.median)<0;
        Stats.median(freq)=rad2deg(stats.median)+360;
    elseif rad2deg(stats.median)>0;
        Stats.median(freq)=rad2deg(stats.median);
    end
    % SKEWNESS 
    Stats.skewness(freq)=stats.skewness0;
    % KURTOSIS
    Stats.kurtosis(freq)=stats.kurtosis0;
    
    
    
    %     Stats.Meansz(freq)=nanmean(degreez+(degreez<0)*360); % AVERAGE DEGREES
    % Working
    % load('WorkSpace2debugFMA.mat') %\/\/\/\/\/\/
    %
    % lfp=[EEG_DownSampledTimestamps',EEG_DownSampledData'];
    %
    % theta = FilterLFP(lfp,'passband','theta');
    %
    % phi = Phase(theta);
    %
    % positions=[(data_video_smoothfilt(:,2)-min(data_video_smoothfilt(:,2)))/range(data_video_smoothfilt(:,2)),...
    % (data_video_smoothfilt(:,3)-min(data_video_smoothfilt(:,3)))/range(data_video_smoothfilt(:,3))];
    %
    % [data,stats] = PhasePrecession(positions,spks_VEL_working(:,1),phi)
    %
 
end
end

% ts_spike_working = interp1(EEG_DownSampledTimestamps', EEG_DownSampledTimestamps', spks_VEL(:,1), 'nearest');
%
% [~,ia,~] = intersect(EEG_DownSampledTimestamps', ts_spike_working);
% N = size(EEG_DownSampledTimestamps');
% spikeVec = zeros(N(1,1),1);
% spikeVec(ia) = 1;
% EEG_DownSampledTimestamps_working = [EEG_DownSampledTimestamps' EEGthetaData' spikeVec(:,1)];
%
% spks_VEL_working = EEG_DownSampledTimestamps_working(EEG_DownSampledTimestamps_working(:,3) == 1,:);
%
%
%
% mints=(spks_VEL_working(1,1));
% maxts=(spks_VEL_working(end,1));
% startts=find(EEG_DownSampledTimestamps_working==mints);
% endts=find(EEG_DownSampledTimestamps_working==maxts);
%
% figure(5)
% plot(EEG_DownSampledTimestamps_working((startts:endts),1)',EEGthetaData(1,(startts:endts)))
% hold on
% scatter(spks_VEL_working(:,1)',spks_VEL_working(:,2)')
%
% figure(6); subplot(2,1,1)
% plot(EEG_DownSampledTimestamps_working((startts:endts),1)',EEGthetaData(1,(startts:endts)))
% hold on
% subplot(2,1,2)
% y=ones(length(EEG_DownSampledTimestamps_working(:,1)));
% stem(spks_VEL_working(:,1)')
%
% ts_spike_working = interp1(EEG_DownSampledTimestamps', EEG_DownSampledTimestamps', spks_VEL(:,1), 'nearest');
%
% ThetaPh_working=rad2deg(ThetaPh);
%
% [~,ia,~] = intersect(EEG_DownSampledTimestamps', ts_spike_working);
% N = size(EEG_DownSampledTimestamps');
% spikeVec = zeros(N(1,1),1);
% spikeVec(ia) = 1;
% EEG_DownSampledTimestamps_working = [EEG_DownSampledTimestamps' ThetaPh_working' spikeVec(:,1)];
%
% spks_VEL_working = EEG_DownSampledTimestamps_working(EEG_DownSampledTimestamps_working(:,3) == 1,:);

% mints=(spks_VEL_working(1,1));
% maxts=(spks_VEL_working(end,1));
% startts=find(EEG_DownSampledTimestamps_working==mints);
% endts=find(EEG_DownSampledTimestamps_working==maxts);

% figure(6)
% plot(EEG_DownSampledTimestamps_working((startts:endts),1)',ThetaPh_working(1,(startts:endts)))
% hold on
% scatter(spks_VEL_working(:,1)',spks_VEL_working(:,2)')

% figure(7)
% Convertedto360=spks_VEL_working(:,2)+(spks_VEL_working(:,2)<0)*360;
% histogram(Convertedto360,length(spks_VEL_working))
% figure(8)
% raster(Convertedto360)
% MeanResultantPhase=circ_r(deg2rad(Convertedto360));
% figure(9)
% rose(Convertedto360)
%
% %
% x = cos(pi/4*(0:99));
% y = hilbert(x);
% sigphase = (unwrap(angle(y)))';
% X = ones(length(sigphase),2);
% X(:,2) = (1:length(sigphase))';
% beta = X\sigphase;
% beta(2)
%
% %%
% x = cos(pi/4*(0:100));
% y = hilbert(x);
% sigphase = atan2(imag(y),real(y));
% X = ones(length(sigphase),2);
% X(:,2) = (1:length(sigphase))';
% beta = X\sigphase;
% beta(2)
% %%
% figure (6)
% plot(ts_seconds,abs(EEGthetaData).^2)

