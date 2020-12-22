% create .avi with neuralynx vt1.nvt timestamp

[Timestamps,Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 0 0 0 1 1], 1, 1, [] );
[x,y,color] = ExtractFromPoints(Targets);

%%
% create the video writer with 1 fps
writerObj = VideoWriter('targets.avi');
% open the video writer
open(writerObj);
for i = 1:10000 %length(Timestamps)
    PlotFrameFromTarget( x,y,color,i)
    writeVideo(writerObj, getframe(gcf))
end

% close the writer object
close(writerObj);
%%


function PlotFrameFromTarget(x,y,color,i)

my_im = zeros(480,720,3);
hold on;   % make sure all targets are plotted in the same window

x_coordinate = x(:,i);    % x coordinate for plotting target
y_coordinate = y(:,i);    % y coordinate for plotting target

y_coordinate = 480 - y_coordinate;  % scale the y value due to a upper left origin for the video tracker.

red_idx = find(color(:,i,1));
green_idx = find(color(:,i,2));
my_im(y_coordinate,x_coordinate,red_idx) = 1;
my_im(y_coordinate,x_coordinate,green_idx) = 1;

imshow(my_im)
end
