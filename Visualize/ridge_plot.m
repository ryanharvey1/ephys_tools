function ridge_plot(varargin)
% ridge_plot: ridge_plot plots are partially overlapping line plots 
% that create the impression of a mountain range. They can be quite 
% useful for visualizing changes in distributions over time or space.
%
% Ryan H 2019

p = inputParser;
p.addParameter('x',peaks); % nxm matrix where each n represents each line
p.addParameter('color_map',0); % if you want a viridis color gradient
p.addParameter('FaceColor',[.3 .3 .3]); % if you want the same color for every line
p.addParameter('FaceAlpha',1); % transparency 
p.addParameter('EdgeColor',[0 0 0]); % transparency 
p.addParameter('LineWidth',.5); % LineWidth
p.parse(varargin{:});

x = p.Results.x;
color_map = p.Results.color_map;
FaceColor = p.Results.FaceColor;
FaceAlpha = p.Results.FaceAlpha;
EdgeColor = p.Results.EdgeColor;
LineWidth = p.Results.LineWidth;

if color_map
    colors=colormap(viridis(length(x)));
end
for i=size(x,1):-1:1
    h = area(x(i,:)+i*2);
    if color_map
        h.FaceColor = colors(i,:);
    else
        h.FaceColor = FaceColor;
    end
    h.FaceAlpha = FaceAlpha;
    h.EdgeColor = EdgeColor;
    h.LineWidth = LineWidth;
    hold on
end
end