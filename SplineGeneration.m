clc
clear all

%% Setting up the spline
init_time = [0,2,7,12]';                 % Initial time vector
q = [0, 2, 12, 5]';                      % Initial points vector
qdot_b = [0, 0]';                        % Spline boundary velocities
qdot_max = 3;                            % Maximum velocity that can be reached
h0 = diff(init_time);                    % Knots vector 
w = abs(q(2:end) - q(1:end-1))/qdot_max; % Lower bound on knots 

%% Setting up the spline knots optimization
options = optimoptions(@fmincon,'Algorithm','sqp');
[h,fval] = fmincon(@objfun, h0, [], [], [], [], w, h0,... 
   @(h)confuneq(h, q, qdot_b, qdot_max), options);

%% Resulting spline based on the optimized knots with given constraints.
[s0, s1, s2, s3] = cubic_spline(h, q, qdot_b); % Generating spline constants
new_time = zeros(length(h)+1, 1); % New knot vector
new_time (1) = 0;
for i=2:length(h)+1
    new_time(i) = new_time (i-1) + h(i-1);
end
plot(new_time, q, 'o');
hold on
plot_cubic_spline(new_time, s0, s1, s2, s3);