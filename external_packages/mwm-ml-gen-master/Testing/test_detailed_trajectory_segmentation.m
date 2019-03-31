
mini = Inf;
maxi = -1;
c = [];
speed = 2;

for i = 1:length(segmentation_configs.TRAJECTORIES.items)

    trajectory_index = i;
    trajectory_length = 99999;

    trajectory_points = segmentation_configs.TRAJECTORIES.items(trajectory_index).points;
    [seg_points,distances,rem_points,rem_distances] = detailed_trajectory_segmentation(trajectory_points, trajectory_length);
    
    a = min(distances(:,1));
    b = max(distances(:,1));
    c = [c;mean(distances(:,1))];
    if a < mini
        mini = a;
    end
    if b > maxi
        maxi = b;
    end

end
disp(a);
disp(b);
disp(min(c));
disp(max(c));

%%

trajectory_index = 1;
trajectory_length = 99999;

trajectory_points = segmentation_configs.TRAJECTORIES.items(trajectory_index).points;
[seg_points,distances,rem_points,rem_distances] = detailed_trajectory_segmentation(trajectory_points, trajectory_length);

f1 = figure;
for i = 1:length(seg_points)
    plot_arena(segmentation_configs);
    hold on;
    plot_trajectory(segmentation_configs.TRAJECTORIES.items(trajectory_index));
    hold on
    points = seg_points{i};
    plot(points(1,2), points(1,3), 'r*', 'LineWidth', 2);
    for j = 2:size(points,1)
        plot(points(j-1:j,2), points(j-1:j,3), '-', 'LineWidth', 2, 'Color', 'red', 'LineStyle', '-');
    end
    plot(points(end,2),points(end,3), 'ro', 'LineWidth', 2);
    k=1;
    pause(speed);
    cla
end
close(f1);