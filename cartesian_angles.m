function [theta, phi] = cartesian_angles(points)
% Convert Cartesian points on the sphere to colatitude and azimuth.

x = points(:, 1);
y = points(:, 2);
z = points(:, 3);

theta = acos(max(-1, min(1, z)));
phi = mod(atan2(y, x), 2 * pi);
end
