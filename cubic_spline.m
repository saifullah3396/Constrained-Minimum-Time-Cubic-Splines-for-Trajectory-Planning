function [s0, s1, s2, s3, mLower, mUpper] = cubic_spline(h, q, qdot_b)

n = length(h) + 1;
d = (q(2:n) - q(1:n-1))./h;

lower = [h(1:n-1);0];
main  = [2*h(1);2*(h(1:n-2) + h(2:n-1));2*h(n-1)];
upper = [0;h(1:n-1)];

T = spdiags([lower main upper], [-1 0 1], n, n);
rhs = [6*(d(1) - qdot_b(1)); 6*(d(2:n-1)-d(1:n-2)); 6*(qdot_b(2) - d(n-1))];
m = T\rhs;
mLower = m(1:end-1);
mUpper = m(2:end);
s0 = mLower./(6.*h);
s1 = mUpper./(6.*h);
s2 = q(2:end)./h - mUpper.*h./(6);
s3 = q(1:end-1)./h - mLower.*h./(6);



