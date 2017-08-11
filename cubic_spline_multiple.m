function [s0, s1, s2, s3, qddotLower, qddotUpper, hrep] = cubic_spline_multiple(nvar, h, q, qdot_b)

%% Dimension check
assert(...
	isequal(size(q), [length(h) + 1, nvar]), ...
    'The trajectory points must be of size (n x nvar)' ...
)
assert(...
    isequal(size(qdot_b), [2, nvar]), ...
    'The boundary velocities must be of size (2 x nvar)' ...
)

%% Setup

n = length(h) + 1;
hrep = repmat(h, [1, nvar]);
d = (q(2:n, :) - q(1:n-1, :)) ./ hrep;

lower = [h(1:n-1); 0];
main  = [2 * h(1); 2 * (h(1:n - 2) + h(2:n - 1)); 2 * h(n - 1)];
upper = [0; h(1:n - 1)];

A = spdiags([lower main upper], [-1 0 1], n, n);
Ainv = inv(A);

for j = 1:nvar
    qddot(:, j) = Ainv * [6*(d(1, j) - qdot_b(1, j)); 6*(d(2:n-1, j)-d(1:n-2, j)); 6*(qdot_b(2, j) - d(n-1, j))];
end
qddotLower = qddot(1:end-1, :);
qddotUpper = qddot(2:end, :);
s0 = qddotLower ./ (6 .* hrep);
s1 = qddotUpper ./ (6 .* hrep);
s2 = q(2:end, :) ./ hrep - qddotUpper.* hrep./(6);
s3 = q(1:end-1, :) ./ hrep - qddotLower.* hrep./(6);



