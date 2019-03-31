function [MI,norm_phaseAmp] = EEG_CFC(data)
%EEG_CFC analyzes neuralynx EEG data for phase coupling between theta and
%low gamma. Uses methods obtained from Tort et al., 2010 and calc_mi from Benjamin Voloh 
%https://github.com/bvoloh/cfc_analysis 

%   INPUT: 
%       - data: data structure compiled from postprocess.m (see function
%       description for details). 
%   OUTPUT: 
%       - MI - modulation index
%       - CV - coefficient of variation (valid if below .10)
%       - norm_phaseAmp 


% Start Analysis ==========================================================

%CROSS FREQUENCY COUPLING ANALYSIS BY TORT ET AL., 2010

%Initialize variables from data structure

EEGthetaData =data.lfp.theta; %Filtered Theta
EEGrawSig=data.lfp.signal; %Raw LFP
Fs=data.lfp.lfpsamplerate; %Sampling rate

[btheta,atheta] = butter(3,[30/(Fs/2) 60/(Fs/2)]); % butterworth filter low gamma normalized by the nyquist frequency

signal_filtered = filtfilt(btheta,atheta,signal);
EEGlgData=data.lfp.signal

EEGlgData = eegfilt(data.lfp.signal,NewSFreq,30,60);  % LOW GAMMA from Tort et al., 2010

%1. GET PHASE ANGLE FOR THETA AND AMPLITUDE FOR LOW GAMMA
LG_amp=abs(hilbert(-EEGlgData)); %instantaneous amplitutde 
ThetaPh=angle(hilbert(-EEGthetaData)); %instantaneous phase

%BIN LOW GAMMA AMPLITUDE INTO 18deg THETA PHASE ANGLES
edges=deg2rad(-180:18:180);
LG_ampIdx=discretize(ThetaPh,edges);

%GET MEAN AMPLITUDE VALUES FOR EACH THETA BIN
phaseAmp=[];
for i=unique(LG_ampIdx)
    phaseAmp(:,i)=nanmean(LG_amp(LG_ampIdx==i));
end

%GET MI USING CODE FROM Benjamin Voloh (based off Tort et al., 2010)
centers=1:20;
x = -3:.01:3;
Q=pdf(makedist('Uniform'),x);
repmat(.1,length(phaseAmp))
MI=calc_mi(centers, phaseAmp, Q);

%NORMALIZE MEAN AMPLITUDE BY DIVIDING EACH BIN VALUE OVER SUM OVER BINS FOR
%PLOTTING
norm_phaseAmp=phaseAmp/sum(phaseAmp);

%GET SHUFFLED DIST 
[modindex_out, powPhsDists_out, bincenters] = computeShuffledModIndex(400,LG_amp,ThetaPh,4)

end

% \\\\\\\\\\\\\\\\\\\\ LOCAL FUNCTIONS BELOW \\\\\\\\\\\\\\\\\\\\

%shuffling from Hasselmo toolbox
function [modindex_out, powPhsDists_out, bincenters] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,thetarange)
  fprintf('Shuffling (%i): ',shuffles);
  modindex = nan(size(gammaamps,1),size(thetaangles,1),shuffles);
  powPhsDists = nan(size(gammaamps,1),36,shuffles);

  for i = 1:shuffles
    % keep user feel confident that all is well
    fprintf('%i ',i);
    
    % find random offset to shift all theta angles by
    offsetOkay = false;
    while ~offsetOkay
      offset = round(rand*size(thetaangles,2));
      offsetOkay = offset>1 && offset<size(thetaangles,2);
    end

    % shift theta angles
    ta = [thetaangles(:,offset+1:end), thetaangles(:,1:offset)];
    % compute mod indices for this shift
    [modindex(:,:,i), powPhsDists(:,:,i), bincenters] = computeModIndex(gammaamps,ta,thetarange);
  end

  % find distro values
  modindex_out(:,:,1) = mean(modindex,3);
  modindex_out(:,:,2) = std(modindex,[],3);
  powPhsDists_out(:,:,1) = mean(powPhsDists,3);
  powPhsDists_out(:,:,2) = std(powPhsDists,[],3);
  fprintf('\n');
end
