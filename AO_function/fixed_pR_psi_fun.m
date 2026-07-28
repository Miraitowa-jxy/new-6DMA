function f_val = fixed_pR_psi_fun(x, v_phaseshift, params)
p_R = x(1:2).';
psi = x(3:5).';
[~,f_val,~,~,~,~,~,~,~,~,~] = bulid_H(p_R, psi, v_phaseshift, params);
if isfield(params, 'phase_objective_scale')
    objective_scale = params.phase_objective_scale;
elseif isfield(params, 'constant_factor')
    objective_scale = params.constant_factor;
else
    objective_scale = 1;
end
f_val = objective_scale * f_val;  % 相位子问题目标缩放
end
