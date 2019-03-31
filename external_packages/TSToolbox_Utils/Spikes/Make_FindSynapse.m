
% USAGE:
% [synTimes,synWeightsZ,synWeightsR,synSig,cch,b,cellPairId] = Make_FindSynapse(S,shank)
% 
% all matrices are post-syn (y) x pre-syn(x)

function [synTimes,synWeightsZ,synWeightsR,synSig,cch,b,cellPairId] = Make_FindSynapse(S,shank)



%Parameters
bin = 0.5;
nbBins = 100;
nbBinsLg = 1000;

cch = [];
cchSame = [];

cellPairId = [];
cellPairIdSame = [];
shameShIx = [];

for x=1:length(S)
    h=waitbar(x/length(S));
    rgx = Range(S{x});
    rx = length(Range(S{x}));
    for y=x+1:length(S)
        
       rgy = Range(S{y});
       [h,b] = CrossCorr(rgx,rgy,bin,nbBins);
       cch = [cch h*rx*bin/1000];
       cellPairId = [cellPairId;[x y]];
     
       if shank(x)==shank(y)
           [h,bl] = CrossCorr(rgx,rgy,bin,nbBinsLg);
           cchSame = [cchSame h*rx*bin/1000];
           cellPairIdSame = [cellPairIdSame;[x y]];
           shameShIx = [shameShIx;1];
       else
           shameShIx = [shameShIx;0];
       end       
    end
end
close(h)

[synLat,synStrZ,synStrR,bounds] = FindSynapse(cch(:,~shameShIx),'synWin',[0 8]);
synTimes = pair2mat(synLat,cellPairId(~shameShIx,:),length(S));
synWeightsZ = pair2mat(synStrZ,cellPairId(~shameShIx,:),length(S));
synWeightsR = pair2mat(synStrR,cellPairId(~shameShIx,:),length(S));
synSig = pair2mat(bounds,cellPairId(~shameShIx,:),length(S));

[synLat,synStrZ,synStrR,bounds] = FindSynapse_SameShank(cchSame,'synWin',[0 8]);
synTimes = pair2mat(synLat,cellPairIdSame,synTimes);
synWeightsZ = pair2mat(synStrZ,cellPairIdSame,synWeightsZ);
synWeightsR = pair2mat(synStrR,cellPairIdSame,synWeightsR);
synSig = pair2mat(bounds,cellPairIdSame,synSig);
