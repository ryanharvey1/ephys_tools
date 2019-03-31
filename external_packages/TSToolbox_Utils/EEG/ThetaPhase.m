function [thetaMod,thetaPh,thetaTrghs,thetaPks,goodEp] = ThetaPhase(fbasename,channel,ep,S,varargin)

% [thetaMod,thetaPh,thetaTrghs,thetaPks] = ThetaPhase(fbasename,channel,ep)
% optional:
% [thetaMod,thetaPh,thetaTrghs,thetaPks] = ThetaPhase(fbasename,channel,ep,folder)
% to indicate another folder where the lfp file is located
% 
% Adrien Peyrache, 2013

highTh = 5000;

if isa(ep,'cell')
    nEp = length(ep);
else
    nEp = 1;
    ept = cell(1,1);
    ept{1} = ep;
    ep = ept;
    clear ept;
end

if ~isempty(varargin)
    folder = varargin{1};
    if ~strcmp(folder(end),'/')
        folder = [folder '/'];
    end
    lfpName = [folder fbasename];
else
    lfpName = fbasename;
end

xml_data = LoadXml([fbasename '.xml']);
lfp = LoadLfp(lfpName,xml_data.nChannels,channel);

dLfp = double(Data(lfp));
dLfp = resample(dLfp,1,5);
rg = Range(lfp);
rg = rg(1:5:end);   

lfp = tsd(rg,abs(dLfp));
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

% Let's look for the troughs. At 250Hz (1250/5), 2 troughs cannont be
% closer than 20 points (if theta reaches 12Hz);
trghs = LocalMinima(dLfp,20,0);
pks = LocalMinima(-dLfp,20,0);
thetaTrghs = tsd(rg(trghs),dLfp(trghs));
thetaPks = tsd(rg(pks),dLfp(pks));

[dumy ix] = Restrict(thetaTrghs,thetaPks,'align','next');
thetaTrghs = subset(thetaTrghs,unique(ix));
[dumy ix] = Restrict(thetaPks,thetaTrghs,'align','prev');
thetaPks = subset(thetaPks,unique(ix));
thetaTrghs = Restrict(thetaTrghs,thetaPks,'align','next');

% Here we interpolate between trough in time. Hilbert is bad (non uniform
% distribution of phases in most cases).
% ph = interp1(Range(thetaTrghs),2*pi*(0:length(Data(thetaTrghs))-1),rg); 
% ph = ph + 0.04*randn(length(ph),1);
% Then, we just need to mod this value with 2pi
% ph = tsd(rg,mod(ph,2*pi));

[w, f, tw, coh, phases] = getWavelet(dLfp,1250/5, 6, 12, 16, 'MORLET','var',0);
[dummy, ix] = max(w);
clear w coh
ix = ix + [0 cumsum(repmat(length(f),[1 size(phases,2)-1]))];
ph = phases(ix);

ph = tsd(tw,ph(:));
thetaPh = cell(length(S),1);
thetaMod = zeros(length(S),3,nEp);

S = Restrict(S,goodEp);

for c=1:length(S)
    thetaPh{c} = Restrict(ph,S{c});

    for e = 1:nEp
        
        ep{e} = intersect(goodEp,ep{e});
        cellPh = Data(Restrict(thetaPh{c},ep{e}));
        if ~isempty(cellPh)
            [mu, Kappa, pval] = CircularMean(cellPh);
            thetaMod(c,:,e) = [mu, pval, Kappa];
        else
            thetaMod(c,:,e) = [NaN NaN NaN];
        end
       
    end
   
    if 0
        figure(1),clf
        hist(cellPh,25)
        keyboard
    end
    
end
        