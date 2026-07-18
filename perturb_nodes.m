function perturbed = perturb_nodes(points, strength)
% Perturb the points tangentially and project them back onto the sphere.

directions = randn(size(points));
directions = directions - sum(directions .* points, 2) .* points;
lengths = vecnorm(directions, 2, 2);
lengths(lengths == 0) = 1;
directions = directions ./ lengths;

perturbed = points + strength * directions;
perturbed = perturbed ./ vecnorm(perturbed, 2, 2);
end
