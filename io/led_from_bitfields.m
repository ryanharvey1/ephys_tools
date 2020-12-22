% create .avi with neuralynx vt1.nvt timestamp

[Timestamps,Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 0 0 0 1 1], 1, 1, [] );
[x,y,color] = ExtractFromBitfields(Points);

%% Write the movie to file
% create the video writer (default 30fps)
writerObj = VideoWriter('VT1.avi');
% open the video writer
open(writerObj);
for i = 1:1000 %length(Timestamps)
   led_frame =  PlotFrameFromTarget( x,y,color,i);
    writeVideo(writerObj, im2frame(led_frame))
end

% close the writer object
close(writerObj);
%%

function [x,y,color] = ExtractFromBitfields( nlx_targets )
[targets,records] = size(nlx_targets);
%initialize the variables
x = zeros(targets,records);
y = zeros(targets,records);
red = zeros(targets,records);
green = zeros(targets,records);

binary_target_string = arrayfun(@(x) dec2bin(x, 32),nlx_targets,'UniformOutput',false);

x = cellfun(@(x) bin2dec(x(21:32)),binary_target_string); % x position
y = cellfun(@(x) bin2dec(x(5:16)),binary_target_string); % y position
red = cellfun(@(x) bin2dec(x(2)),binary_target_string); % pure red
green = cellfun(@(x) bin2dec(x(3)),binary_target_string); % pure green

color = cat(3,red,green);
disp('All records processed.');
end


function my_im = PlotFrameFromTarget(x,y,color,i)
my_im = zeros(481,721,3);
x_coordinate = x(:,i) + 1;    % x coordinate for plotting target add one due to zero based indexing 
y_coordinate = y(:,i) + 1;    % y coordinate for plotting target

temp_color = color(:,i,1:2);
for ii = 1:length(x_coordinate)
my_im(y_coordinate(ii),x_coordinate(ii),1) = temp_color(ii,1,1);
my_im(y_coordinate(ii),x_coordinate(ii),2) = temp_color(ii,1,2);
end
end
