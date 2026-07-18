function [points, theta, phi, weights] = fibonacci_nodes(number_of_points)
% Fibonacci points on the sphere with equal weights.

index = (0:number_of_points - 1)';
z = 1 - 2 * (index + 0.5) / number_of_points;
golden_ratio = (1 + sqrt(5)) / 2;
phi = mod(2 * pi * index / golden_ratio, 2 * pi);
radius = sqrt(max(0, 1 - z.^2));

points = [radius .* cos(phi), radius .* sin(phi), z];
theta = acos(z);
weights = repmat(4 * pi / number_of_points, number_of_points, 1);
end
