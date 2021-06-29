% create .avi with neuralynx vt1.nvt timestamp

[Timestamps,Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 0 0 0 1 1], 1, 1, [] );
[x,y,color] = ExtractFromBitfields(Points);
%[x,y,color,valid_targets] = ExtractFromBitfields(Points);

%% Write the movie to file
% create the video writer (default 30fps)
writerObj = VideoWriter('VT1.avi');
% open the video writer
open(writerObj);
for i = 1:1000 %length(Timestamps)
   led_frame =  PlotFrameFromTarget( x,y,color,valid_targets,record_index ); %( x,y,color,i);
    writeVideo(writerObj, im2frame(led_frame))
end

% close the writer object
close(writerObj);
%%

function [x,y,color] = ExtractFromBitfields(nlx_targets)
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

% function [x,y,color,valid_targets] = ExtractFromTargets( nlx_targets )
% 
%     % find how many recs and how many targets per record
%     [max_targets,max_records] = size( nlx_targets );
% 
%     %output number of rec in file
%     str = sprintf('There are %d records in this file',max_records);
%     disp(str);
% 	
% 	%initialize the variables in case no targets are found
% 	x(1) = 0;
% 	y(1) = 0;
% 	color(1) = 0;
% 	valid_targets(1) = 0;
%     % loop to extract all targets from each record
%     for rec_index = 1:max_records
% 
%         %get record
%         current_record = nlx_targets(:,rec_index);
% 
%         % loop to extract all targets from current record
%         for target_index = 1:max_targets
%             
%             % if the bitfield is equal to zero, then we know there is no valid data for that
%             % field, or the rest of the bitfields in the record.
%             if current_record( target_index ) == 0
%                 break;
%             end
% 
%             % extract the x and y positions and store them to be returned to matlab
%             [x( rec_index, target_index ),y( rec_index, target_index )] = ExtractXY( current_record( target_index ) );
%             [color( rec_index, target_index, 1:7 )] = ExtractColor( current_record( target_index ) );
%             
%         end  %end inner loop
% 
%         % record the number of targets within each rec that contain data
%         valid_targets( rec_index ) = target_index - 1;
%             
% 
%     end  %end loop
%     
%     disp('All records processed.');
% end

%----------------------------------------------------------------------------------------------------------------------
%   This function extracts the x and y coordinates from the bitfield for a given target.  
%----------------------------------------------------------------------------------------------------------------------
function [x, y] = ExtractXY(target_value)  

	binary_target_string = dec2bin(target_value, 32);
	x = bin2dec(binary_target_string(21:32));
	y = bin2dec(binary_target_string(5:16));
  end
  
%----------------------------------------------------------------------------------------------------------------------
%	Extracts color information from a target
%----------------------------------------------------------------------------------------------------------------------
function [color] = ExtractColor(target_value)

	binary_target_string = dec2bin(target_value, 32);
	color(1) = bin2dec(binary_target_string(2)); %pure red
	color(2) = bin2dec(binary_target_string(3)); %pure green
	color(3) = bin2dec(binary_target_string(4)); %pure blue
	color(4) = bin2dec(binary_target_string(18)); %raw red
	color(5) = bin2dec(binary_target_string(19)); %raw green
	color(6) = bin2dec(binary_target_string(20)); %raw blue
	color(7) = bin2dec(binary_target_string(17)); %intensity

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

function frame = PlotFrameFromTarget( x,y,color,valid_targets,record_index )

    num_targets = valid_targets(record_index);  % get the number of targets for the given record that are valid
    
    figure(record_index);   % index the given figure window with the same index for the specified frame
    hold on;   % make sure all targets are plotted in the same window
    
    % loop through number of targets for the given rec
    for target_index = 1:num_targets
    
        x_coordinate = x(record_index,target_index);    % x coordinate for plotting target
        y_coordinate = y(record_index,target_index);    % y coordinate for plotting target
        
        y_coordinate = 480 - y_coordinate;  % scale the y value due to a upper left origin for the video tracker.
        
        str_color = GetColorString(color(record_index,target_index,1:7));   % gets the string of the color to plot the target
     
        plot( x_coordinate, y_coordinate, str_color);   % plot target
 
    end
    
    axis([0 720 0 480]);    % set axis for graph
    hold off;   % release the plot
    frame = getframe(gcf);
end
%----------------------------------------------------------------------------------------------------------------------
%  This function returns a string containg color information for plotting the frame.  Color is chosen in a specific
%   order based on preference.  The color yellow is used to represent raw green due to the absence of a light green color string.
%----------------------------------------------------------------------------------------------------------------------
function [color_string] = GetColorString( color_array )

    if ( color_array(2) ~= 0 )      % pure red
        color_string = 'ro';
    elseif ( color_array(4) ~= 0 )  % pure green
        color_string = 'go';        
    elseif ( color_array(6) ~= 0 )  % pure blue
        color_string = 'bo';
    elseif ( color_array(5) ~= 0 )  % raw blue
        color_string = 'co';
    elseif ( color_array(1) ~= 0 )  % raw red
        color_string = 'mo';
    elseif ( color_array(3) ~= 0 )  % raw green
        color_string = 'yo';
    else                            % black for luminance or default
        color_string = 'ko';
    end
end
