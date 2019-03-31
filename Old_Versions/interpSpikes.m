function [B] = interpSpikes(rawTS,filtmat,SpikeTS )
%interpSpikes interpolates spike times with a TS by position matrix
% 
% Input:    rawTS: Unfiltered contiguous timestamps 
%           filtmat: Filtered non-contiguous timestamps / xy position / headangle etc...
%           SpikeTS: Raw timestamps associated with spikes
% 
% Output:   B: TS by Position by binary spike matrix *
%
%   *EXAMPLE OUTPUT
%
%   TIMESTAMP   X           Y          ANGLE      TIMESTAMP SPIKEBINARY
%   424347346   581.4363    40.2363    222.7093   0         0
%   424382259   580.8858    39.6889    222.1315   0.0333    0
%   424548421   580.3895    39.1132    221.5405   0.06667   1
%
% Ryan E Harvey, Laura Berkowitz MAY 2017
% 
% INTERPOLATE SPIKES TO RAW TS
ts_spike = interp1(rawTS(:,1), rawTS(:,1), SpikeTS, 'nearest');
% FIND MATCHES TO FILTERED MATRIX
ts_spike2=ts_spike(ismember(ts_spike,filtmat(:,1)));
% BUILD MATRIX FROM INTERP SPIKE TIMESTAMPS 
ts_spike2working=[];
for i=1:length(ts_spike2)
    ts_spike2working=[ts_spike2working;filtmat(all(repmat(ts_spike2(i),size(filtmat(:,1),1),1)==filtmat(:,1),2),:)];
end
% ADD ONES 
ts_spike2=[ts_spike2working,ones(size(ts_spike2working,1),1)];
% DELETE TIME STAMPS USED IN ABOVE MATRIX
filtmat(ismember(filtmat(:,1),ts_spike2),:)=[];
% COMBINE TIMESTAMP AND FILTERED MATRIX, THEN SORT 
B=sortrows([ts_spike2;[filtmat,zeros(length(filtmat),1)]],1);

