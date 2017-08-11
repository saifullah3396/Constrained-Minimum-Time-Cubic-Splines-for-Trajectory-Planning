function cubic_driver(num_points)

runge = @(t) t.^2

t = [-2,-1,0,1,2]
q = runge(t)

qdot_b = [0, 0]

[s0,s1,s2,s3] = cubic_spline(t',q', qdot_b');

plot_points = 1000;
% tt = linspace(-7,7,plot_points);
% qq = runge(tt);
plot(t,q, 'o');
hold on
% plot(tt,qq,'g');
hold on;
plot_cubic_spline(t,s0,s1,s2,s3);
% You can see it in action bq running the following at the Matlab prompt
