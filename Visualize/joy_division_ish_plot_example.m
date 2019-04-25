% joy_division_ish_plot_example
% Ryan H
fig=figure;
x=peaks(49);
colors=colormap(viridis(length(x)));
for i=1:length(x)
plot(1:length(x),x(i,:)'+i*2,'LineWidth',2,'Color',colors(i,:))
hold on
end
darkBackground(fig,[0.2 0.2 0.2],[0.7 0.7 0.7])
box off
axis off

export_fig('C:\Users\ryanh\Dropbox\yourfigure2.pdf');

