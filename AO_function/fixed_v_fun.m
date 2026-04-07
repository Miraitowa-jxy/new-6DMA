function F= fixed_v_fun(x, v_phaseshift, params)
p_R = x(1:2).';
psi = x(3:5).';
[~,f_val,~,~,~,~,~,~,~,~,~] = bulid_H(p_R, psi, v_phaseshift, params);
H = params.H;
p_B = params.p_B;
p_U = params.p_U;
p_T = params.p_T;
constant_factor = params.constant_factor;
p_R_real = [p_R;H];

Q = rotation_matrix(psi);
L = Q(:,3);
% 计算惩罚项
penalty_lambda = 10; %惩罚系数
value1 = -L.' * (p_R_real - p_B);
value2 = -L.' * (p_R_real - p_U);
value3 = -L.' * (p_R_real - p_T);
penalty = max(0, value1)^2 + max(0, value2)^2 + max(0, value3)^2;
% 计算适应性函数 F
F = f_val + penalty_lambda * penalty;
F = constant_factor * F; % 乘以常数constant_factor,防止matlab因为精度问题导致错误
end

