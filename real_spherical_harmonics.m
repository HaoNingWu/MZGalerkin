function values = real_spherical_harmonics(n, theta, phi)
% Real spherical harmonics of degree at most n.
% Columns are ordered by ell = 0,...,n and m = -ell,...,ell.

theta = theta(:);
phi = phi(:);
number_of_points = numel(theta);
values = zeros(number_of_points, (n + 1)^2);
x = cos(theta);

column = 1;
for ell = 0:n
    associated_legendre = legendre(ell, x');

    for order = -ell:ell
        absolute_order = abs(order);
        polynomial = associated_legendre(absolute_order + 1, :)';
        normalization = sqrt((2 * ell + 1) / (4 * pi) * ...
            factorial(ell - absolute_order) / factorial(ell + absolute_order));

        if order == 0
            values(:, column) = normalization * polynomial;
        elseif order > 0
            values(:, column) = sqrt(2) * normalization * polynomial .* ...
                cos(absolute_order * phi);
        else
            values(:, column) = sqrt(2) * normalization * polynomial .* ...
                sin(absolute_order * phi);
        end
        column = column + 1;
    end
end
end
