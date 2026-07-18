%PLOT_SAMPLING_POINTS Reproduce the two point sets shown in Fig. 1.

number_of_points = 500;
perturbation = 0.2;

rng(1);
points = fibonacci_nodes(number_of_points);
perturbed_points = perturb_nodes(points, perturbation);

figure_handle = figure('Color', 'w', 'Name', 'Sampling points');
layout = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
draw_points(points);
title('Fibonacci');

nexttile;
draw_points(perturbed_points);
title('perturbed');

title(layout, sprintf('Point sets on the sphere, m = %d', number_of_points));

function draw_points(points)
[x, y, z] = sphere(80);
surf(x, y, z, 'FaceColor', [0.72 0.82 0.72], 'FaceAlpha', 0.75, ...
    'EdgeColor', 'none');
hold on;
scatter3(points(:, 1), points(:, 2), points(:, 3), 9, ...
    [0.65 0.05 0.05], 'filled');
axis equal off;
view(135, 24);
camlight headlight;
lighting gouraud;
end
