function tort = computeMI_hilbert(self,varargin)
% Computes the phase modulation of gamma by theta
% using the method described by Tort et al (2010) J Neurophys.
% In short, it computes the
% entropy in the histogram of powers over phases.
%
%
% INPUT ARGS
%  root - cmb object with active_lfp set and lfp loaded
%  filtType - default = [], otherwise 'theta' or 'vel' to process data only
%             with high theta power or fast running speeds.
%  filtParams - default = [], lower and upper bounds to use in the
%               filtering of the data as mentioned above
%  shuffles - default 0, use if you want to compute null dist of modulation
%  thetaRange - default = [4:0.25:12] frequencies to look at phases of
%  gammaRange - defualt = [30:1:120] frequencies to look at the power of
%  ifPlot - Default = 0. If 1, plots in new figure. If anything else, then
%           assumes input is an axis to plot onto
%
%
% OUTPUT ARGS
%  modindex - matrix of size [nG x nT] where nT is the number of theta frequencies
%             analyzed and nG is the number of gamma frequencies analysed
%  thetarange - same as inputed, usable as labels for the matrix when plotting
%  gammarange - same as inputed, usable as labels for the matrix when plotting
%
% [modindex, thetarange, gammarange, powPhsDists, bincenters] = thetaModGamma(root,filtType,filtParams,shuffles,thetarange,gammarange)


p = inputParser;
p.addParamValue('filtType', []);
p.addParamValue('filtParams',   []);
p.addParamValue('shuffles', 0);
p.addParamValue('thetarange',   4:0.25:12);
p.addParamValue('gammarange',   30:1:120);
p.addParamValue('nbins',   36);
p.addParamValue('ifPlot',   0);

p.parse(varargin{:});

thetarange = p.Results.thetarange;
gammarange = p.Results.gammarange;
filtType = p.Results.filtType;
filtParams = p.Results.filtParams;
shuffles = p.Results.shuffles;
nbins = p.Results.nbins;
ifPlot = p.Results.ifPlot;

% get gamma amplitudes
gamma_filtered = GammaFilter(self,gammarange);
gammaamps_raw=abs(hilbert(gamma_filtered));

% get theta filtered signals and phase
thetaangles_raw = extractThetaPhase(self,'hilbert',thetarange);

% select out high quality data (no saturations and no flat lines + apply requested filtering
[gammaamps,thetaangles] = cleanData(self,gammaamps_raw,thetaangles_raw,thetarange,filtType,filtParams);

% double check that there are at least 200 cycles of theta (7Hz), the min
% mentioned by Tort et al (2010).
if size(thetaangles,2)<200*(self.lfpsamplerate/8)
    warning('There may be too few theta cycles!');
end

% compute the values of interest
[tort.modindex, tort.MeanAmp] = computeModIndex(gammaamps,thetaangles,nbins);

% compute significance thresholds for each frequency pairing
[tort.MIthres, tort.modindex_shuffled, tort.MeanAmp_shuffled] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,nbins);

if ifPlot
    
    fig=figure;
    fig.Color=[1 1 1];
    subplot(1,2,1)
    bar([tort.MeanAmp./sum(tort.MeanAmp) tort.MeanAmp./sum(tort.MeanAmp)])
    
    title(['Data MI  (maxmod = ' num2str(max(tort.modindex(:)))])
    xlabel('Theta phase (10 degree bins)'), ylabel('Normalized Gamma Amps')
    
    subplot(1,2,2)
    bar([mean(tort.MeanAmp_shuffled,3)./sum(mean(tort.MeanAmp_shuffled,3))...
         mean(tort.MeanAmp_shuffled,3)./sum(mean(tort.MeanAmp_shuffled,3))])
    
    title(['Surrogate MI (MI threshold 95% cutoff = ' num2str(max(tort.MIthres(:)))])
    xlabel('Theta phase (10 degree bins)'), ylabel('Normalized Gamma Amps')
    
end
end
%###################### LOCAL FUNCTIONS ###################################

%% Cleandata
function [gammaamps,thetaangles] = cleanData(self,gammaamps,thetaangles,thetarange,filtType,filtParams)

% compute specific filters if needed
if exist('filtType','var') && ~isempty(filtType)
    
    % find good amps
    switch(filtType)
        case 'theta'
            thetaamps = abs(hilbert(buttfilt(self.signal,minmax(thetarange),self.lfpsamplerate,'bandpass',4)));
            excInds = thetaamps < max(thetaamps)*filtParams(1) | thetaamps > max(thetaamps)*filtParams(2);
        case 'vel'
            hdvel = interp1(self.tsVel,self.vel,self.ts);
            excInds = hdvel < filtParams(1) | hdvel > filtParams(2);
        otherwise
            error('Unknown filtType');
    end
    excInds = reshape(excInds,1,[]);
else
    excInds = zeros(1,size(gammaamps,2));
end

%%  drop saturated or flatlined data

maxSig = max(self.signal);
minSig = min(self.signal);
t0 = [0-1e-5 0+1e-5]; % thresholds for no dV/dt

[~,bad1]=ThresholdBandDetect_v2(diff(self.signal),t0(1),t0(2),1,5); % flatline data
[~,bad2]= OverThresholdDetect_v2(self.signal,maxSig,1,5); % railed out at top
[~,bad3]= UnderThresholdDetect_v2(self.signal,minSig,1,5); % railed out at bottom

badInds =sum([bad1, bad2, bad3],2)>0; % railed out at bottom

% perform filter
gammaamps(:,excInds|badInds') = [];
thetaangles(:,excInds|badInds') = [];
end

% computes MI with hilbert method for estimation of amplitude and phase via
% tort 2010
function [modindex, meanamps] = computeModIndex(gammaamps,thetaangles,nbins)
% bin gamma amplitudes by theta phases

%Bin gamma amplitude into phase bins of width=binwidth (default 17 degrees)
meanamps=Findmeanamps(thetaangles,gammaamps,nbins);

%compute shannon entropy
shannonent=-sum((meanamps/sum(meanamps)).*log((meanamps/sum(meanamps))));

%modindex is equal to the KL distance of the observed amplitude
%distribution from the uniform distribution divided by log(N).
modindex = (log(nbins)-shannonent) / log(nbins);
end

% %%% computes surrogate analysis for with hilbert method for estimation of amplitude and phase via
% tort 2010 (surrogate data uses phase shifts of theta phase).
function [MIthreshold, modindex_out, MeanAmp_out] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,nbins)
% fprintf('Shuffling (%i): ',shuffles);
rng shuffle
offset=randi((size(thetaangles,2)-2),1,200)+1;

for i = 1:shuffles
    % keep user feel confident that all is well
%      shift theta angles
    ta = [thetaangles(:,offset(1,i)+1:end), thetaangles(:,1:offset(1,i))];
    
    % compute mod indices for this shift
    [modindex, MeanAmp] = computeModIndex(gammaamps,ta,nbins);
    MeanAmp_out(:,:,i)=MeanAmp;
    modindex_out(1,i)=modindex;
end

% find distro values
MIthreshold=prctile(modindex_out,95);
MeanAmp_out(:,:,1) = mean(MeanAmp_out,3);
% fprintf('\n');
end

%
function [meanamps]=Findmeanamps(slowphase,fastamp,nbins)
%this function finds the normalized mean amplitude of the fast wave for each phase of the
%slow wave
%it recieves the phase of the slow wave (in radians from 0-2pi) and the
%corresponding amplitude of the fast wave (in mV)

%it returns a matrix with the normalized mean amplitude (row 2) corresponding to the
%phases(row 1-- max of corresponding phase), and the middle of the phase
%for a bar graph(row 3)

bin=360/nbins;

%BIN LOW GAMMA AMPLITUDE INTO 10deg THETA PHASE ANGLES
edges=deg2rad(-180:bin:180);
LG_ampIdx=discretize(slowphase,edges);

%GET MEAN AMPLITUDE VALUES FOR EACH THETA BIN
meanamps=[];
for i=unique(LG_ampIdx)
    meanamps(:,i)=nanmean(fastamp(LG_ampIdx==i));
end
% add in bin labsls

end

%%
function [phase,pow]=multiphasevec(f,S,Fs,width)
% FUNCTION [phase,pow]=multiphasevec(f,S,Fs,width)
%
% Returns the phase and power as a function of time for a range of
% frequencies (f).
%
% Simply calls phasevec in a loop.
%
% INPUT ARGS:
%   f = [2 4 8];   % Frequencies of interest
%   S = dat;       % Signal to process
%   Fs = 256;      % Sampling frequency
%   width = 6;     % Width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   phase- Phase data [freqs,time]
%   power- Power data [freqs,time]
%

pow = nan(length(f),length(S));
phase = nan(length(f),length(S));

for a=1:length(f)
    [phase(a,:),pow(a,:)]=phasevec(f(a),S,Fs,width);
end

end

%%
function thetaPhs = extractThetaPhase(self,method,band)

if isempty(self)
    error('Please load lfp with active_lfp before calling');
end

switch(lower(method))
    
    case 'waveform'
        if exist('band','var') || ~isempty(band)
            warning('band not used when computing phase by waveform')
        end
        
        %  method described by
        broadBand = [1 60];
        filtSig = buttfilt(self.signal,broadBand,self.lfpsamplerate,'bandpass',4);
        
        % define window detection params
        min_sep = self.lfpsamplerate / 20; % 50ms separation between peak and next trough
        min_length = min_sep; % same for the separation between trough and next peak
        
        % Find troughs and peaks
        % as onsets and offsets of signal with positive slope
        trphsAndPeaks = OverThresholdDetect_v2(sign(diff(filtSig)),0,min_sep,min_length);
        
        % Find zero crossings as sign changes in signal
        upAndDownZeroXings = OverThresholdDetect_v2(sign(filtSig),0,min_sep,min_length);
        
        waveLandmarks = sort([trphsAndPeaks, upAndDownZeroXings]);
        
        error('This code is not complete')
        % this method was to be based off of the methods described by Belluscio
        % et al 2012, however they had a large number of channels
        % simultaneously recorded from which they computed the median signal
        % and then computed the wavefrom from this.  This likely dramatically
        % reduces the number of local minima and maxima found and thereby makes
        % this a practical endeavor.  To make this work on a signal channel of
        % data, the signal would likely have to be filtered more thereby
        % reducing the value of using this method.
        
    case 'wavelet'
        % returns a vector of values, spanning theta range, per time point
        if exist('band','var') & ~isempty(band)
            thetarange = band;
        else
            thetarange = 4:0.1:12;
        end
        [thetaPhs,~] = multiphasevec(thetarange,self.signal,self.lfpsamplerate,16);
        
        
    case 'hilbert'
        % returns a single phase estimate by time point using bandpass filtering
        if exist('band','var') & ~isempty(band)
            thetarange = band;
        else
            thetarange = [4 12];
        end
        
        signal_filtered = ThetaFilter(self,thetarange);
        
        
        thetaPhs = angle(hilbert(signal_filtered));
        
    otherwise
        error('Unknown method')
end
end

%%
function y = morlet(f,t,width)
% function y = morlet(f,t,width)
%
% Morlet's wavelet for frequency f and time t.
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet.
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY
%
% Ole Jensen, August 1998

sf = f/width;
st = 1/(2*pi*sf);
A = 1/sqrt(st*sqrt(pi));
y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);
end

%%
function [y, pow] = phasevec(f,s,Fs,width)
% FUNCTION [y,pow] = phasevec(f,s,Fs,width)
%
% Return the phase as a function of time for frequency f.
% The phase is calculated using Morlet's wavelets.
%
% INPUT ARGS:
%   f = 8;     % Frequency of interest
%   s = dat;   % The eeg signal to evaluate
%   Fs = 256;  % Sample rate
%   width = 6; % width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   y- Phase as a function of time
%   pow- Power as a function of time
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);

y = conv(s,m);

pow = abs(y).^2;
pow = pow(ceil(length(m)/2):length(pow)-floor(length(m)/2));

l = find(abs(y) == 0);
y(l) = 1;

% normalizes phase estimates to length one
y = y./abs(y);
y(l) = 0;

y = angle( y(ceil(length(m)/2):length(y)-floor(length(m)/2)) );
end

function signal_filtered = ThetaFilter(self,thetarange)
% Takes 'signal' and bandpasses it to theta frequencies (4 to 12 Hz)
%
% Arguments
%
% signal - arbitrary signal to be filtered
% Fs - sample frequency
%
% signal_filtered = ThetaFilter(signal, Fs)

Wn_theta = [thetarange(1)/(self.lfpsamplerate/2) thetarange(end)/(self.lfpsamplerate/2)]; % normalized by the nyquist frequency

[btheta,atheta] = butter(3,Wn_theta);

signal_filtered = filtfilt(btheta,atheta,self.signal);
end

function signal_filtered = GammaFilter(self,gammarange)
% Takes 'signal' and bandpasses it to gamma frequency
%
% Arguments
%
% self - structure with self.signal (arbitrary signal to be filtered) and
% self.lfpsamplerate (sample frequency)
% gammarange - range of gamma frequency of interest
%


Wn_gamma = [gammarange(1)/(self.lfpsamplerate/2) gammarange(end)/(self.lfpsamplerate/2)]; % normalized by the nyquist frequency

[bgamma,agamma] = butter(3,Wn_gamma);

signal_filtered = filtfilt(bgamma,agamma,self.signal);
end


%% CODE GRAVEYARD


%%
% function [MI,MeanAmp] = ModIndex_v1(thetaangles,gammaamps,position,fig)
%
% % the eegfilt routine employed below is obtained from the EEGLAB toolbox
% % (Delorme and Makeig J Neurosci Methods 2004)
%
%
% Phase=thetaangles; % this is getting the phase time series
%
% Amp=gammaamps; % getting the amplitude envelope
%
% % Now we search for a Phase-Amp relation between these frequencies by
% % caclulating the mean amplitude of the AmpFreq in each phase bin of the
% % PhaseFreq
%
% % Computing the mean amplitude in each phase:
%
% nbin=length(position);
% winsize = 2*pi/nbin;
%
% MeanAmp=zeros(1,nbin);
% for j=1:nbin
% I = find(Phase <  position(j)+winsize & Phase >=  position(j));
% MeanAmp(j)=mean(Amp(I));
% end
% % the center of each bin (for plotting purposes) is position+winsize/2
%
% % quantifying the amount of amp modulation by means of a
% % normalized entropy index (Tort et al PNAS 2008):
%
% MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);
%
% end