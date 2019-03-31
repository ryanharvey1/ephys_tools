function TS=temporalstability(spkts,fullts)
% temporalstability: check to see if spike events reoccur over time
% Important for filtering out potential cells that might just be head bump or
% unplug noise not otherwise filtered by cluster quality measures. Such
% data can seem very spatially modulated when it is truely not. 
% 
% Input:
%           spkts: timestamps for ever spike
%           fullts: every timestamp 
% 
% Output:   TS: 0 to 1
%               1: stable over time
%               0: only a single spiking event...probably noise
% 
% Ryan Harvey 2018


% bin spike times to assess stability (1 bin = 10 second)
edges=linspace(min(fullts),max(fullts),round(length(fullts)/300));
temporalhist=histcounts(spkts,edges);

TS=sum(temporalhist>0)/length(temporalhist);


% % check to see if spiking event occurs over time
% if sum(temporalhist>max(temporalhist)*.20)>2
%     TS=1;
% else
%     TS=0;
% end
end
