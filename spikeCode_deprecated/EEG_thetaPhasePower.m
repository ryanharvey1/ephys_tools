% EEG_thetaPhase filters raw Neuyralynx CSC files and calculates a filtered
% theta signal and phase
% 
% Input: 
%       - eegfile = path to CSC file
% 
% Output: 
%       - EEGreshaped = raw EEG; 
%       - EEGtsUpSampled = raw upsampled timestamps 
%       - EEGthetaData = filtered by theta EEG
%       - ThetaPh = theta phase
% 
% created by Ben C; patched together from Nlx2MatCSC from Neuralynx and McN lab code CalcThetaPhase.m, EEGtsUpSample.m,
% InterpDownSample.m, RestrictEEG.m
% 
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\spikeCode'));  
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\Cell analysis'));  
addpath(genpath('F:\Users\BClarkLab\Documents\MATLAB\RyanDirectory\chronux_2_11')); 
% load eeg file
eegfile = 'F:\Users\reharvey\Place_Cell_Data\PAE_Rat\LS17\2017-03-02_12-13-31\CSC2.ncs';

% extract CSC data with mex
[Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples]= Nlx2MatCSC(eegfile, [1 1 1 1 1], 0, 1); 

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

% quick check that downsampling works
% LengthDS = length(EEG_DownSampledTimestamps);
% dsSampleTS = EEG_DownSampledTimestamps(round(LengthDS/2):round(LengthDS/2)+NewSFreq-1);
% dsSampleData = EEG_DownSampledData(round(LengthDS)/2:round(LengthDS/2)+NewSFreq-1);
% SampleInd = find(dsSampleTS(1)<=EEGtsUpSampled & EEGtsUpSampled<=dsSampleTS(end));
% SampleTS = EEGtsUpSampled(SampleInd);
% SampleData = EEGreshaped(SampleInd);
% figure
% plot(SampleTS,SampleData, 'k')
% hold on
% plot(dsSampleTS, dsSampleData, 'r')

% filter between 4 and 12 Hz.
EEGthetaData = eegfilt(EEG_DownSampledData,NewSFreq,4,12);

% determine theta phase using the Hilbert transform
ThetaPh = angle(hilbert(-EEGthetaData));

% create x- axis (time in msec)
ts_seconds = index/1000;

% plot raw EEG
figure (1);
subplot(5,1,1), plot(ts_seconds,EEG_DownSampledData, 'k');
axis tight
title('Raw EEG');
ylabel('Voltage (mV)')
xlabel('Time (sec)')
% xlim([0 60]);
hold on
% plot theta filtered EEG
subplot(5,1,2), plot(ts_seconds,EEGthetaData, 'r');
axis tight
title('Theta Filtered EEG (6-10 Hz)');
ylabel('Voltage (mV)')
xlabel('Time (sec)')
% xlim([0 60]);

% Compute and Plot Power Spectrum, code from Rhythms of the Cortex SfN 2009 Short Course #2
v0 = EEG_DownSampledData .* hann(length(EEG_DownSampledData))'; %Multiply data by Hann taper
pow = (abs(fft(v0)).^2) / length(v0); %Compute the spectrum.
pow = 10*log10(pow/max(pow)); %Convert to decibels.
pow = pow(1:length(ts_seconds)/2+1); %Ignore negative frequencies.
dt = 0.001; %Define the time resolution.
f0 = 1/dt; %Determine the sampling frequency.
df = 1/max(ts_seconds); %Determine the frequency resolution.
fNQ = f0/2; %Determine the Nyquist frequency.
faxis = (0:df:fNQ); %Construct frequency axis.
subplot(5,1,3), plot(faxis, pow);
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

% Create Periodogram
% dt = 0.001; %The sampling interval in seconds
Fs = 1/dt; %The sampling frequency in Hertz  
subplot(5,1,4), periodogram(EEG_DownSampledData,[],length(EEG_DownSampledData),Fs); %function from Matlab signal processing toolbox  
xlim([0 20]);

% Compute and plot power spectrum, function from Chronux toolbox
% dt = 0.001; %The sampling interval in seconds
% Fs = 1/dt; %The sampling frequency in Hertz                                   
params.Fs = Fs;                            
params.tapers = [1 2]; %Chronux parameters                    
[S,f] = mtspectrumc(EEG_DownSampledData, params); %function from Chronux toolbox  
S = 10*log10(S/max(S));
subplot(5,1,5), plot(f,S);
xlim([0 20])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

% % Plot spectrogram
% % Window size is 1000 pts = 1 s, overlap is 500 pts = 0.5 s, compute fft over 1000 pts
% [S,F,T] = spectrogram(EEG_DownSampledData,1000,500,1000,Fs); 
% S = abs(S);
% figure (2)
% subplot(2,1,1), imagesc(T,F,10*log10(S/max(S(:)))); colorbar;
% axis xy
% ylim([0 20])
% xlabel('Time (s)')
% ylabel('Freq (Hz)')
% hold on
% 
% % Compute the multi-taper spectrogram using Chronux toolbox
% params.Fs = Fs;                       
% params.tapers = [2 1]; %Chronux specific parameters                             
% [Smtm,Tmtm,Fmtm] = mtspecgramc(EEG_DownSampledData,[1 0.5],params);    
% subplot(2,1,2), imagesc(Tmtm, Fmtm, 10*log10(Smtm/max(Smtm(:)))'); colorbar;
% axis xy
% ylim([0 20])
% xlabel('Time (s)')
% ylabel('Freq (Hz)')