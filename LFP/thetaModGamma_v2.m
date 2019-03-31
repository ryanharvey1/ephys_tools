function CFC= thetaModGamma_v2(self,varargin)
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
p.addParamValue('thetarange',   6:0.25:10);
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


% self.signal=self.signal(active_lfp_ind,:);

% initialize outputs

if ~exist('partial','var')||isempty(partial)
    % get gamma amplitudes
    [~,gammaamps] = multiphasevec(gammarange,self.signal,self.lfpsamplerate,8);
    
    % get theta filtered signals and phase
    thetaangles = extractThetaPhase(self,'wavelet',thetarange);
end

% select out high quality data (no saturations and no flat lines + apply requested filtering
[gammaamps,thetaangles] = cleanData(self,gammaamps,thetaangles,thetarange,filtType,filtParams);

% double check that there are at least 200 cycles of theta (7Hz), the min
% mentioned by Tort et al (2010).
if size(thetaangles,2)<200*(self.lfpsamplerate/8)
    warning('There may be too few theta cycles!');
end


% COMPUTE MI COMUDULATION
[CFC.modindex, CFC.MeanAmp, ~] = computeModIndex(gammaamps,thetaangles,thetarange,nbins);
CFC.maxMI=max(CFC.modindex(:));

if shuffles>0
    % COMPUTE SIGNIFICANCE THRESHOLDS FOR EACH FREQUENCY PAIRING
    [CFC.modindex_shuffled, CFC.modindex_shufSTD] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,thetarange,nbins);
    CFC.maxMI_shuffled=max(CFC.modindex_shuffled(:));
else
    CFC.modindex_shuffled=NaN; CFC.modindex_shufSTD=NaN;
    CFC.maxMI_shuffled=NaN;
end


% STORE INPUT VARS
CFC.thetarange=thetarange;
CFC.gammarange=gammarange;
CFC.nbins=nbins;
CFC.fs=self.lfpsamplerate;


if ifPlot
    if ifPlot==1
        figure;
        ax = gca;
    else
        ax = ifPlot;
    end
    axes(ax);
%     subplot(1,2,1)
    imagesc(thetarange,gammarange,CFC.modindex);
    axis xy
    
    title(['Data MI  (maxmod = ' num2str(max(CFC.modindex(:)))])
    xlabel('Theta Frequency'), ylabel('Gamma Frequency')
    
%     subplot(1,2,2)
%     imagesc(thetarange,gammarange,CFC.modindex_shuffled);
%     axis xy
%     
%     title(['Surrogate MI  (maxmod = ' num2str(max(CFC.modindex_shuffled(:)))])
%     xlabel('Theta Frequency'), ylabel('Gamma Frequency')
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

%
function [modindex, meanamps, bincenters] = computeModIndex(gammaamps,thetaangles,thetarange,nbins)
% bin gamma amplitudes by theta phases
meanamps = nan(nbins,size(gammaamps,1),size(thetaangles,1));
bincenters = nan(nbins,size(thetaangles,1));
for i=1:length(thetarange)
    [meanamps(:,:,i),bincenters]=Findmeanamps(thetaangles(i,:)',gammaamps',nbins);
end

% get back to the way things used to be
meanamps = permute(meanamps,[2 3 1]);

% normalize the area to one
meanampsNorm = meanamps ./ repmat(sum(meanamps,3),[1,1,size(meanamps,3)]);

% compute the mod indexes
shannonent = sum(meanampsNorm .* log10(meanampsNorm),3); %Shannon entropy of a dist.

%compute shannon entropy
modindex = (log10(nbins)+shannonent) ./ log10(nbins);

% ELN edit 120512 - The following code will change powPhsDists to be the
% power distrobution over theta phases for the theta band with the
% maximum modulation index.

%%% compute mean phase alignments over window of high modulation
% powPhsDists = meanamps;


end

% %%
function [modindex_shuffled, modindex_std] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,thetarange,nbins)
% fprintf('Shuffling (%i): ',shuffles);
modindex = nan(size(gammaamps,1),size(thetaangles,1),shuffles);
powPhsDists = nan(size(gammaamps,1),nbins,shuffles);

% find random offset to shift all theta angles by
rng shuffle
offset=randperm((size(thetaangles,2)-2),200)+1;

parfor i = 1:shuffles
    % keep user feel confident that all is well
%     fprintf('%i ',i);
    
    % shift theta angles
    ta = [thetaangles(:,offset(1,i)+1:end), thetaangles(:,1:offset(1,i))];
    
    % compute mod indices for this shift
    [modindex_out(:,:,i), ~, ~] = computeModIndex(gammaamps,ta,thetarange,nbins);
    
end

% find distro values
modindex_shuffled= nanmean(modindex_out,3);
modindex_std= nanstd(modindex_out,[],3);
%   powPhsDists_out(:,:,1) = mean(powPhsDists,3);
%   powPhsDists_out(:,:,2) = std(powPhsDists,[],3);
% fprintf('\n');

end

%
function [meanamps, bincenters]=Findmeanamps(slowphase,fastamp,nbins)
%this function finds the normalized mean amplitude of the fast wave for each phase of the
%slow wave
%it recieves the phase of the slow wave (in radians from 0-2pi) and the
%corresponding amplitude of the fast wave (in mV)

%it returns a matrix with the normalized mean amplitude (row 2) corresponding to the
%phases(row 1-- max of corresponding phase), and the middle of the phase
%for a bar graph(row 3)



%%make into column vectors if not
if size(slowphase,1)==1
    slowphase=slowphase';
end

% double check data matches
if size(fastamp,1)~=size(slowphase,1)
    error('fastamp must have time in rows');
end

% define phase bins
sorting = linspace(-pi, pi, nbins+1);

% perform sorting
meanamps = nan(length(sorting)-1,size(fastamp,2));
for bin = 1:length(sorting)-1
    currInds = slowphase > sorting(bin) & slowphase <= sorting(bin+1);
    meanamps(bin,:) = mean(fastamp(currInds,:),1);
end

% add in bin labsls
bincenters = mean([sorting(1:end-1);sorting(2:end)])'+pi;

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
        thetaPhs = angle(hilbert(buttfilt(self.signal,thetarange,self.lfpsamplerate,'bandpass',4)));
        
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