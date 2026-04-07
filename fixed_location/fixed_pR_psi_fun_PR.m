function f_val = fixed_pR_psi_fun_PR(x, v_phaseshift, params)
p_R = params.p_R.';
psi = x.';
[~,f_val,~,~,~,~,~,~,~,~,~] = bulid_H(p_R, psi, v_phaseshift, params);
constant_factor = params.constant_factor;
f_val = constant_factor * f_val;  % 乘以常数constant_factor,防止matlab因为精度问题导致错误
end
