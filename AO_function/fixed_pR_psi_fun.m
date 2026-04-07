function f_val = fixed_pR_psi_fun(x, v_phaseshift, params)
p_R = x(1:2).';
psi = x(3:5).';
[~,f_val,~,~,~,~,~,~,~,~,~] = bulid_H(p_R, psi, v_phaseshift, params);
constant_factor = params.constant_factor;
f_val = constant_factor * f_val;  % 乘以常数constant_factor,防止matlab因为精度问题导致错误
end
