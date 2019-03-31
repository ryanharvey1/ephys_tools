function raster(X,samplerate)
% raster: simple raster plot
% Input:
%           X: cell array of spike timestamps
%           sample rate: sample rate of video (to scale x axis)
%
% example 
% figure;raster(data.Spikes,data.samplerate)
%
% Ryan E Harvey 2018

[~,I] = sort(cellfun(@length,X));
X = X(I);


for i=1:length(X)
   plot(X{i},zeros(length(X{i}),1)+i,'.k') ;hold on
end
ylim([0 length(X)])
xlabel('Time(sec)')
ylabel('Cells')
ax=gca;
set(ax,'XTickLabel',round(ax.XTick/samplerate))
box off
end