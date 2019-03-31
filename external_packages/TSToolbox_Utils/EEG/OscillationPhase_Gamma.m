function [gammaMod,gammaPh,goodEp] = GammaPhase(lfp,S,ep,freqBand)

% [gammaMod,gammaPh] = GammaPhase(fbasename,channel,ep,freqBand)
% INPUTS:
%   channel  : LFP channel number (base 1)
%   ep       : intervalSet object identifying where the signal should be looked at
%   freqBand : 2 value vector indicating low and high frequancy range of gamma band
% optional:
%  [gammaMod,gammaPh] = GammaPhase(fbasename,channel,ep,freqBandfolder)
%  to indicate another folder where the lfp file is located
% 
% Adrien Peyrache, 2013

highTh      = 5000;
lowNoiseFq  = 200; %cut-off frequency to identify noise epochs in LFP

if numel(freqBand) ~=2
    error('input argument ''freqBand'' must have 2 values')
elseif freqBand(1)>=freqBand(2)    
    error('input argument ''freqBand'' should have increasing values')
end

if isa(ep,'cell')
    nEp = length(ep);
else
    nEp = 1;
    ept = cell(1,1);
    ept{1} = ep;
    ep = ept;
    clear ept;
end

%Downsampling lfp for faster computation
dLfp    = double(Data(lfp));
dLfp    = resample(dLfp,1,2);
rg      = Range(lfp);
rg      = rg(1:2:end);
lfp     = tsd(rg,abs(dLfp));

%Removing artifacts
badEp = thresholdIntervals(lfp,highTh); %empirical values here to remove artifacts
badEp = mergeCloseIntervals(badEp,2);
badEp = intervalSet(Start(badEp)-0.9999,End(badEp)+0.9999);
badEp = mergeCloseIntervals(badEp,2);
badEp = dropShortIntervals(badEp,0.001);
goodEp = timeSpan(lfp)-badEp;
goodEp = dropShortIntervals(goodEp,1);
clear lfp

%Amp. limit for pks/trgh detection
% lfp = tsd(rg,dLfp);
% ampL = std(Data(Restrict(lfp,goodEp)));

%NoiseDetection
noisePow    = gaussianFilter(dLfp,1250/2,[lowNoiseFq 1250/2]);
noisePow    = gaussFilt(abs(noisePow),4*1250/lowNoiseFq);
noiseMed    = median(noisePow);
noiseMedVar = median(abs(noisePow-noiseMed));
noiseMedZ   = abs((noisePow-noiseMed)/noiseMedVar);

cleanEp     = thresholdIntervals(tsd(rg,noiseMedZ(:)),7,'Direction','Below');

goodEp      = intersect(goodEp,cleanEp);

[w, f, tw, coh, phases] = getWavelet(dLfp,1250/2, freqBand(1), freqBand(2), 16);
[maxVal, ix] = max(w);
powerBand = mean(w,2);

clear w coh
ix = ix + [0 cumsum(repmat(length(f),[1 size(phases,2)-1]))];
ph = phases(ix);

ph = tsd(tw,ph(:));

gammaPh = cell(length(S),1);
gammaMod = zeros(length(S),3,nEp);

S = Restrict(S,goodEp);

for c=1:length(S)
    gammaPh{c} = Restrict(ph,S{c});

    for e = 1:nEp
        
        ep{e} = intersect(goodEp,ep{e});
        cellPh = Data(Restrict(gammaPh{c},ep{e}));
        if ~isempty(cellPh)
            [mu, Kappa, pval] = CircularMean(cellPh);
            gammaMod(c,:,e) = [mu, pval, Kappa];
        else
            gammaMod(c,:,e) = [NaN NaN NaN];
        end
       
    end
   
    if 0
        figure(1),clf
        hist(cellPh,25)
        keyboard
    end
    
end
        