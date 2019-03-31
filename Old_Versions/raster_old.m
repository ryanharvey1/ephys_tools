function raster(spikes)

for i=1:length(spikes)
    in=spikes{i}';
    
    plot([in;in],[ones(size(in));zeros(size(in))]+(i-1),'k-');hold on
    
end
set(gca,'TickDir','out') % draw the tick marks on the outside
set(gca,'YTick', []) % don't draw y-axis ticks
% set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
set(gca,'Color',get(gcf,'Color')) % match figure background
set(gca,'YColor',get(gcf,'Color')) % hide the y axis
box off
end

% rasterplot

spikes=S;

for i=1:length(spikes)
lengths(i)=length(spikes{i});
end
[s,I]=sort(lengths);

figure;
ii=1;
for i=I
   plot(spikes{i},zeros(length(spikes{i}),1)+ii,'.k');hold on 
   ii=ii+1;
end
ylabel('Cells')
xlabel('time')