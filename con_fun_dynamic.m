function [qdot_con1, qdot_con2] = con_fun_dynamic(h, n_discretized, ndof, q, qdot_b, qdot_max, dyn_params_M, dyn_params_B, dyn_params_C, dyn_params_G)

t = zeros(length(h)+1, 1);
t(1) = 0;
for i=2:length(h)+1
t(i) = t(i-1) + h(i-1);
end
tUpper = repmat(t(2:end), [1, ndof]);
tLower = repmat(t(1:end-1), [1, ndof]);
qdot_cont1 = zeros(length(h) * ndof, 1);
vC2 = zeros(length(h) * ndof, 1);

[s0, s1, s2, s3, qddotLower, qddotUpper, hrep] = cubic_spline_multiple(ndof, h, q, qdot_b); 
qdot_endpoints = -3 .* s0' *(hrep).^2 + s2' - s3';
qdot_con1_temp = abs(qdot_endpoints) - repmat(qdot_max, [1, ndof]);
qdot_con1 = reshape(qdot_con1_temp, [numel(qdot_con1_temp), 1]);
inflexion_check = (qddotLower .* qddotUpper) < 0;
t_mid = tLower + hrep .* qddotLower ./ (qddotLower - qddotUpper)
qdot_midpoint = -3 .* s0' * (tUpper - t_mid).^2 + 3 .* s1' *(t_mid - tLower).^2 + ... 
                s2' - s3';
qdot_con2_temp = inflexion_check .* (abs(qdot_midpoint) - repmat(qdot_max, [1, ndof]));
qdot_con2 = reshape(qdot_con2_temp, [numel(qdot_con2_temp), 1]);

% tauC = zeros((n_discretized * length(h) + 1) * ndof, 1);
% t_index = 1;
% for m = 1:ndof
%     d_index = 1;
%     for i=1:length(h)
%         diff_time = h(i)/n_discretized;
%         for t_inc = tLower(i):diff_time:(tUpper(i)-diff_time)
%             qddot = -6*s0(i)*(tUpper(i) - t_inc) + 6*s1(i)*(t_inc-tLower(i)) + ... 
%             tauC(t_index, 1) = dyn_params_M{d_index} * qddot;
%             t_index = t_index + 1;
%             d_index = d_index + 1;
%         end
%     end
%     tauC(t_index, 1) = t_index;
%     t_index = t_index + 1;
% end



