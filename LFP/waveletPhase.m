% waveletPhase tries to replicate methods used in Amemiya & Redish (2018)
% and Belluscio et al (2012). 
%
%Inputs: 
%       Dependent on data structure used in B.Clark lab (2017-2018)
%Script procedure: 
%       1. Loads lfp for subject/parses data   
%       2. Broad bandpass filter 
%       3. Filter for theta
%       4. Find local minima (0°) /maxima (180°), Ascending points
%       (90°) & descending points (270°)
%       5. Loads lfp for subject 


%parse inputs 
EEGrawSig=data.lfp.signal; %Raw LFP
Fs=data.lfp.lfpsamplerate; %Sampling rate


for i=1:length(EEGrawSig,2)
%Bandpass filter(1–80 Hz)
bpFilt = BandpassFilter(EEGrawSig, Fs, [1 80]);

% Theta Filter 
thetaFilt = ThetaFilter(bpFilt, Fs);

%Find local minima/maxima for waveform based phase extraction & theta
%asymmetry - Theta waveform asymmetry was quantified by the “asymmetry index,”
% log[length of ascending part of theta wave] / log[length of descending part of theta wave]

%local minima and maxima in the theta range (6–12 Hz) were determined as peaks and trough of theta, respectively

end