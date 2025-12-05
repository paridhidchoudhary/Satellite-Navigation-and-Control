function plot_multiple_results(time, all_actual, all_actual_del, vapf_all, initial_conditions)
    colors = lines(length(all_actual)); % distinguishable colors

    % 3D Trajectories APF
    figure;
    hold on;
    for i = 1:length(all_actual)
        plot3(all_actual{i}(4,:), all_actual{i}(5,:), all_actual{i}(6,:), 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    scatter3(goal(1), goal(2), goal(3), 100, 'r', 'filled');
    title('APF Trajectories - Varying Initial Conditions');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legendStrings = arrayfun(@(i) sprintf('Start #%d', i), 1:length(all_actual), 'UniformOutput', false);
    legend([legendStrings, {'Goal'}]);
    grid on; axis equal;

    % Similarly, add for D-APF
    figure;
    hold on;
    for i = 1:length(all_actual_del)
        plot3(all_actual_del{i}(4,:), all_actual_del{i}(5,:), all_actual_del{i}(6,:), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    scatter3(goal(1), goal(2), goal(3), 100, 'r', 'filled');
    title('D-APF Trajectories - Varying Initial Conditions');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend([legendStrings, {'Goal'}]);
    grid on; axis equal;
end
