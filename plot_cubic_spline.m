function plot_cubic_spline(t,s0,s1,s2,s3)

n = length(t);

inner_points = 100;

for i=1:n-1
   tt = linspace(t(i),t(i+1),inner_points);
   ti = repmat(t(i),1,inner_points);
   ti1 = repmat(t(i+1),1,inner_points);
   qq = s0(i)*(ti1 - tt).^3 + s1(i)*(tt-ti).^3 + ... 
     s2(i)*(tt-ti) + s3(i)*(ti1 - tt);
   qqdot = -3*s0(i)*(ti1 - tt).^2 + 3*s1(i)*(tt-ti).^2 + ... 
     s2(i) - s3(i);
   if(tt == ti)
     qqdot
   end
   plot(tt,qq,'b')
   plot(tt,qqdot, 'r');
%    plot(t(i),0,'r');
end
% Here is a function that constructs a cubic spline and plots in on the famous Runge function:
