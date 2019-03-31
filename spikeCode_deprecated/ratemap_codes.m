%% plot filled and smoothed rate map 
NaNsForSmoothedR = NaN(1,nBins+1);
NaNsForSmoothedC = NaN(nBins,1);
SmoothRateMapGraph = [SmoothRateMap,NaNsForSmoothedC;NaNsForSmoothedR];
figure(2), subplot(1,1,1), h = pcolor(SmoothRateMapGraph);
axis off square tight
hold on
box off
set(h, 'EdgeColor', 'none');

%% plot spatial autocorr
NaNsForAutocorrR = NaN(1,nBins*2);
NaNsForAutocorrC = NaN(nBins*2-1,1);
SpatialAutoCorrGraph = [SpatialAutoCorr,NaNsForAutocorrC;NaNsForAutocorrR];
figure(3), subplot(1,1,1), h = pcolor(SpatialAutoCorrGraph);
axis off square tight
hold on
box off
set(h, 'EdgeColor', 'none');

[9:09] 
better code for the rate maps