function [degree, order] = basis_indices(n)
% Indices of real spherical harmonics through degree n.

number_of_functions = (n + 1)^2;
degree = zeros(number_of_functions, 1);
order = zeros(number_of_functions, 1);

column = 1;
for ell = 0:n
    block = column:(column + 2 * ell);
    degree(block) = ell;
    order(block) = (-ell:ell)';
    column = block(end) + 1;
end
end
