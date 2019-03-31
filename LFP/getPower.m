function power= getPower(self,varargin)
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
p.addParamValue('thetarange',   4:0.25:12);
p.addParamValue('gammarange',   25:1:100); %broad range


p.parse(varargin{:});

thetarange = p.Results.thetarange;
gammarange = p.Results.gammarange;
filtType = p.Results.filtType;
filtParams = p.Results.filtParams;

% get gamma amplitudes
% [~,gammaamps] = multiphasevec(gammarange,self.signal,self.lfpsamplerate,6);

%Get amps via tort 2010
gamma_filtered = GammaFilter(self,gammarange);
gammaamps=abs(hilbert(gamma_filtered));

% get theta filtered signals and phase
[thetaangles,self.thetaPow] = extractThetaPhase(self,'hilbert',thetarange);

% select out high quality data (no saturations and no flat lines + apply requested filtering
[gammaamps,thetaangles,self] = cleanData(self,gammaamps,thetaangles,thetarange,filtType,filtParams);

% power.maxgamma=max(gammaamps(:));

%Define the track quadrants 
minLim=self.xlimits(1,1); maxLim=self.xlimits(1,2);
x=rescale([self.xcoord;minLim;maxLim],0,1);
x(end-1:end)=[];

new_x = interp1(self.tsVel,x,self.ts);

quad1IDX= new_x >= 0 & new_x <= .25;
quad2IDX= new_x >= .26 & new_x <= .5;
quad3IDX= new_x >= .51 & new_x <= .75;
quad4IDX= new_x >= .76 & new_x <= 1;

nbins=36;
shuffles=200;

%Whole lap
power.thetaPowLap=self.thetaPow(:,quad1IDX');
power.gammaPowLap=gammaamps(:,quad1IDX');
[power.MIlap, power.MIlap_amps] = computeModIndex(gammaamps,thetaangles,nbins);
[power.MIlap_shuffled, power.MIlap_out, power.MIlap_amps_shuffled]=computeShuffledModIndex(shuffles,gammaamps,thetaangles,nbins); %hilbert method


% %QUADRANT 1
power.thetaPowq1=self.thetaPow(:,quad1IDX');
power.gammaPowq1=gammaamps(:,quad1IDX');
[power.MIq1, power.MIq1_amps] = computeModIndex(gammaamps(:,quad1IDX'),thetaangles(:,quad1IDX'),nbins);
[power.MIq1_thres_shuffled, power.MIq1_out, power.MIq1_amps_shuffled]=computeShuffledModIndex(shuffles,gammaamps(:,quad1IDX'),thetaangles(:,quad1IDX'),nbins); %hilbert method

% [MIq1_shuffled]=computeShuffledModIndex(shuffles,gammaamps(:,quad1IDX'),...
%thetaangles(:,quad1IDX'),thetarange,nbins); %morlet method - MI over
%individual frequencies
% logi=MIq1_shuffled<power.MIq1;
% power.MIq1pval=1-((sum(logi,3)/200));
% power.meanMIq1_shuffled=mean(MIq1_shuffled,3);

%QUADRANT 2
power.thetaPowq2=self.thetaPow(:,quad2IDX');
power.gammaPowq2=gammaamps(:,quad2IDX');
[power.MIq2, power.MIq2_amps] = computeModIndex(gammaamps(:,quad2IDX'),thetaangles(:,quad2IDX'),nbins);
[power.MIq2_thres_shuffled, power.MIq2_out, power.MIq2_amps_shuffled]=computeShuffledModIndex(shuffles,gammaamps(:,quad2IDX'),thetaangles(:,quad2IDX'),nbins);

% [MIq2_shuffled]=computeShuffledModIndex(shuffles,gammaamps(:,quad2IDX'),thetaangles(:,quad2IDX'),thetarange,nbins);
% logi=MIq2_shuffled<power.MIq2;
% power.MIq2pval=1-((sum(logi,3)/200));
% power.meanMIq2_shuffled=mean(MIq2_shuffled,3);

%QUADRANT 3
power.thetaPowq3=self.thetaPow(:,quad3IDX');
power.gammaPowq3=gammaamps(:,quad3IDX');
[power.MIq3, power.MIq3_amps] = computeModIndex(gammaamps(:,quad3IDX'),thetaangles(:,quad3IDX'),nbins);
[power.MIq3_thres_shuffled, power.MIq3_out, power.MIq3_amps_shuffled]=computeShuffledModIndex(shuffles,gammaamps(:,quad3IDX'),thetaangles(:,quad3IDX'),nbins);

% [MIq3_shuffled]=computeShuffledModIndex(shuffles,gammaamps(:,quad3IDX'),thetaangles(:,quad3IDX'),thetarange,nbins);
% logi=MIq3_shuffled<power.MIq3;
% power.MIq3pval=1-((sum(logi,3)/200));
% power.meanMIq3_shuffled=mean(MIq3_shuffled,3);

%QUADRANT 4
power.thetaPowq4=self.thetaPow(:,quad4IDX');
power.gammaPowq4=gammaamps(:,quad4IDX');
[power.MIq4, power.MIq4_amps] = computeModIndex(gammaamps(:,quad4IDX'),thetaangles(:,quad4IDX'),nbins);
[power.MIq4_thres_shuffled, power.MIq4_out, power.MIq4_amps_shuffled]=computeShuffledModIndex(shuffles,gammaamps(:,quad4IDX'),thetaangles(:,quad4IDX'),nbins);

% [MIq4_shuffled]=computeShuffledModIndex(shuffles,gammaamps(:,quad4IDX'),thetaangles(:,quad4IDX'),thetarange,nbins);
% logi=MIq4_shuffled<power.MIq4;
% power.MIq4pval=1-((sum(logi,3)/200));
% power.meanMIq4_shuffled=mean(MIq4_shuffled,3);

end

%###################### LOCAL FUNCTIONS ###################################

%% Cleandata
function [gammaamps,thetaangles,self] = cleanData(self,gammaamps,thetaangles,thetarange,filtType,filtParams)

% compute specific filters if needed
if exist('filtType','var') && ~isempty(filtType)
    
    % find good amps
    switch(filtType)
        case 'theta'
            theta_filtered = ThetaFilter(self,thetarange);
            thetaamps = abs(hilbert(theta_filtered));
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
self.ts(:,excInds|badInds') = [];
self.thetaPow(:,excInds|badInds') = [];
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

if isempty(gammaamps) || length(gammaamps)<143
    MIthreshold=NaN;
    modindex_out=NaN;
    MeanAmp_out=NaN;
    return
end
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
% function [modindex, meanamps, bincenters] = computeModIndex(gammaamps,thetaangles,thetarange,nbins)
% % bin gamma amplitudes by theta phases
% meanamps = nan(nbins,size(gammaamps,1),size(thetaangles,1));
% bincenters = nan(nbins,size(thetaangles,1));
% for i=1:length(thetarange)
%     [meanamps(:,:,i),bincenters]=Findmeanamps(thetaangles(i,:)',gammaamps',nbins);
% end
% 
% % get back to the way things used to be
% meanamps = permute(meanamps,[2 3 1]);
% 
% % normalize the area to one
% meanampsNorm = meanamps ./ repmat(sum(meanamps,3),[1,1,size(meanamps,3)]);
% 
% % compute the mod indexes
% shannonent = sum(meanampsNorm .* log10(meanampsNorm),3); %Shannon entropy of a dist.
% 
% %compute shannon entropy
% modindex = (log10(nbins)+shannonent) ./ log10(nbins);
% 
% % ELN edit 120512 - The following code will change powPhsDists to be the
% % power distrobution over theta phases for the theta band with the
% % maximum modulation index.
% 
% %%% compute mean phase alignments over window of high modulation
% % powPhsDists = meanamps;
% 
% 
% end
% 
% % %%
% function [modindex_shuffled, modindex_std] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,thetarange,nbins)
% modindex = nan(size(gammaamps,1),size(thetaangles,1),shuffles);
% powPhsDists = nan(size(gammaamps,1),nbins,shuffles);
% 
% % find random offset to shift all theta angles by
% rng shuffle
% offset=randi((size(thetaangles,2)-2),1,200)+1;
% 
% parfor i = 1:shuffles
%     % keep user feel confident that all is well
% %     fprintf('%i ',i);
%     
%     % shift theta angles
%     ta = [thetaangles(:,offset(1,i)+1:end), thetaangles(:,1:offset(1,i))];
%     
%     % compute mod indices for this shift
%     [modindex_out(:,:,i), ~, ~] = computeModIndex(gammaamps,ta,thetarange,nbins);
%     
% end
% 
% % find distro values
% modindex_shuffled= modindex_out;
% % modindex_std= std(modindex_out,[],3);
% %   powPhsDists_out(:,:,1) = mean(powPhsDists,3);
% %   powPhsDists_out(:,:,2) = std(powPhsDists,[],3);
% % fprintf('\n');
% 
% end

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
%%
function [thetaPhs,thetaPow] = extractThetaPhase(self,method,band)

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
        [thetaPhs,thetaPow] = multiphasevec(thetarange,self.signal,self.lfpsamplerate,16);
        
        
    case 'hilbert'
        % returns a single phase estimate by time point using bandpass filtering
        if exist('band','var') & ~isempty(band)
            thetarange = band;
        else
            thetarange = [4 12];
        end
        
        signal_filtered = ThetaFilter(self,thetarange);
        
        
        thetaPhs = angle(hilbert(signal_filtered));
        thetaPow = abs(hilbert(signal_filtered));

        
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
% 
% function [meanamps, bincenters]=Findmeanamps(slowphase,fastamp,nbins)
% %this function finds the normalized mean amplitude of the fast wave for each phase of the
% %slow wave
% %it recieves the phase of the slow wave (in radians from 0-2pi) and the
% %corresponding amplitude of the fast wave (in mV)
% 
% %it returns a matrix with the normalized mean amplitude (row 2) corresponding to the
% %phases(row 1-- max of corresponding phase), and the middle of the phase
% %for a bar graph(row 3)
% 
% 
% 
% %%make into column vectors if not
% if size(slowphase,1)==1
%     slowphase=slowphase';
% end
% 
% % double check data matches
% if size(fastamp,1)~=size(slowphase,1)
%     error('fastamp must have time in rows');
% end
% 
% % define phase bins
% sorting = linspace(-pi, pi, nbins+1);
% 
% % perform sorting
% meanamps = nan(length(sorting)-1,size(fastamp,2));
% for bin = 1:length(sorting)-1
%     currInds = slowphase > sorting(bin) & slowphase <= sorting(bin+1);
%     meanamps(bin,:) = mean(fastamp(currInds,:),1);
% end
% 
% % add in bin labsls
% bincenters = mean([sorting(1:end-1);sorting(2:end)])'+pi;
% 
% end

% %% CODE GRAVEYARD


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