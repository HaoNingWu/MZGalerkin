%RUN_EXPERIMENT Reproduce the spherical Galerkin experiment from the paper.
n = 8;
gamma = 0.3;
point_counts = [100 120 180 260 380 560 820 1200 1800 2600 ...
    3600 5000 7500 10000 15000]';
perturbation = 0.08;

rng(1);

[degree, ~] = basis_indices(n);
dimension = numel(degree);
operator_eigenvalues = 1 + gamma * degree .* (degree + 1);

exact_coefficients = zeros(dimension, 1);
number_of_modes = min(12, dimension);
exact_coefficients(1:number_of_modes) = randn(number_of_modes, 1);
exact_coefficients = exact_coefficients / norm(exact_coefficients);
rhs = operator_eigenvalues .* exact_coefficients;

fprintf('Fibonacci points\n');
clean = run_family(point_counts, n, operator_eigenvalues, rhs, ...
    exact_coefficients, 0);

fprintf('\nPerturbed Fibonacci points\n');
perturbed = run_family(point_counts, n, operator_eigenvalues, rhs, ...
    exact_coefficients, perturbation);

disp(clean);
disp(perturbed);

results.clean = clean;
results.perturbed = perturbed;
results.degree = n;
results.gamma = gamma;
results.perturbation = perturbation;
results.figure = plot_results(clean, perturbed);


function values = run_family(point_counts, n, eigenvalues, rhs, exact, strength)
number_of_cases = numel(point_counts);
error = zeros(number_of_cases, 1);
gram_defect = zeros(number_of_cases, 1);
gram_condition = zeros(number_of_cases, 1);

for k = 1:number_of_cases
    number_of_points = point_counts(k);
    [points, theta, phi, weights] = fibonacci_nodes(number_of_points);

    if strength > 0
        points = perturb_nodes(points, strength);
        [theta, phi] = cartesian_angles(points);
    end

    harmonics = real_spherical_harmonics(n, theta, phi);
    gram = harmonics' * (weights .* harmonics);
    galerkin_matrix = gram * diag(eigenvalues);
    coefficients = galerkin_matrix \ rhs;

    error(k) = norm(coefficients - exact);
    gram_defect(k) = norm(gram - eye(size(gram)), 2);
    gram_condition(k) = cond(gram);

    fprintf(['m = %5d   error = %.3e   ||G-I||_2 = %.3e   ' ...
        'cond(G) = %.3e\n'], number_of_points, error(k), ...
        gram_defect(k), gram_condition(k));
end

values = table(point_counts, error, gram_defect, gram_condition, ...
    'VariableNames', {'m', 'error', 'gram_defect', 'gram_condition'});
end

function figure_handle = plot_results(clean, perturbed)
clean_color = [0.30 0.40 0.55];
perturbed_color = [0.70 0.40 0.25];

figure_handle = figure('Color', 'w', 'Name', 'MZ Galerkin experiment');
layout = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
draw_pair(clean.m, clean.error, perturbed.m, perturbed.error, ...
    clean_color, perturbed_color);
xlabel('number of points $m$', 'Interpreter', 'latex');
ylabel('$L^2$ error $\|u_n^{(m)}-u_n\|_{L^2}$', 'Interpreter', 'latex');
title('(a) Galerkin error');

nexttile;
draw_pair(clean.m, clean.gram_defect, perturbed.m, perturbed.gram_defect, ...
    clean_color, perturbed_color);
xlabel('number of points $m$', 'Interpreter', 'latex');
ylabel('Gram defect $\|G_m-I\|_2$', 'Interpreter', 'latex');
title('(b) Inner-product defect');

nexttile;
draw_pair(clean.gram_defect, clean.error, perturbed.gram_defect, ...
    perturbed.error, clean_color, perturbed_color);
xlabel('Gram defect $\|G_m-I\|_2$', 'Interpreter', 'latex');
ylabel('$L^2$ error $\|u_n^{(m)}-u_n\|_{L^2}$', 'Interpreter', 'latex');
title('(c) Error versus defect');
legend('Fibonacci', 'perturbed', 'Location', 'northwest');

title(layout, 'Quadrature-based Galerkin approximation on the sphere');
end

function draw_pair(x1, y1, x2, y2, color1, color2)
loglog(x1, y1, 'o-', 'Color', color1, 'LineWidth', 1.4, ...
    'MarkerSize', 5, 'MarkerFaceColor', 0.9 * color1);
hold on;
loglog(x2, y2, 's--', 'Color', color2, 'LineWidth', 1.4, ...
    'MarkerSize', 5, 'MarkerFaceColor', 0.9 * color2);
grid on;
box on;
set(gca, 'FontSize', 11);
end
