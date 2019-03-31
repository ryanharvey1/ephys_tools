function [cellMod,cellPh] = OscillationPhase_Ripples(lfp,S,ep)

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

dLfp = double(Data(lfp));
rg = Range(lfp);

[w, f, tw, coh, phases] = getWavelet(dLfp,1250, 100, 200, 16);
[dummy, ix] = max(w);
clear w coh
ix = ix + [0 cumsum(repmat(length(f),[1 size(phases,2)-1]))];
ph = phases(ix);

ph = tsd(tw,ph(:));
cellPh = cell(length(S),1);
cellMod = zeros(length(S),3,nEp);

%S = Restrict(S,goodEp);

for c=1:length(S)
    cellPh{c} = Restrict(ph,S{c});

    for e = 1:nEp
        
        tmpPh = Data(Restrict(cellPh{c},ep{e}));
        if ~isempty(tmpPh)
            [mu, Kappa, pval] = CircularMean(tmpPh);
            cellMod(c,:,e) = [mu, pval, Kappa];
        else
            cellMod(c,:,e) = [NaN NaN NaN];
        end
       
    end
   
    if 0
        figure(1),clf
        hist(cellPh,25)
        keyboard
    end
    
end
        