function [vC1, vC2] = confuneq(h, q, qdot_b, qdot_max)

t = zeros(length(h)+1, 1);
t(1) = 0;
for i=2:length(h)+1
t(i) = t(i-1) + h(i-1);
end

[s0, s1, s2, s3, mLower, mUpper] = cubic_spline(h, q, qdot_b);

tUpper = t(2:end);  
tLower = t(1:end-1);

vC1 = zeros(length(h), 1);
vC2 = zeros(length(h), 1);
for i=1:length(h)
    qdot = -3*s0(i)*(h(i)).^2 + s2(i) - s3(i);
    vC1(i,:) = abs(qdot) - qdot_max;
    if(mLower(i)*mUpper(i) < 0);
        t_star = t(i) + h(i)*mLower(i)/(mLower(i) - mUpper(i))
        qdot_star = -3*s0(i)*(tUpper(i) - t_star).^2 + 3*s1(i)*(t_star-tLower(i)).^2 + ... 
            s2(i) - s3(i);
        vC2(i,:) = abs(qdot_star) - qdot_max;
    end
end