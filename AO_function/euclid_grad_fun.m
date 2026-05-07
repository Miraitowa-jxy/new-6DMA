% function euclid_grad = euclid_grad_fun(p_R, psi, v_phaseshift, params)
% [~,~, hc_H, ht_H, hr, Uc, a_bar_BU_H, Ut, a_bar_BT_H, Ur, a_bar_TB] = bulid_H(p_R, psi, v_phaseshift, params);
% %将目标函数分为三个部分
% f1 = hc_H * ht_H';
% f2 = -hr' * hr;
% f3 = ht_H * hc_H';
% %对各部分函数求导
% delta_f1 = conj(Ut) * a_bar_BU_H.' + conj(Ut) * Uc.' * v_phaseshift;  
% delta_f2 = -Ur' * a_bar_TB - Ur' * Ur * v_phaseshift;
% delta_f3 = conj(Uc) * a_bar_BT_H.' + conj(Uc) * Ut.' * v_phaseshift;
% euclid_grad = f2*f3*delta_f1 + f3*f1*delta_f2 + f1*f2*delta_f3;
% constant_factor = params.constant_factor;
% euclid_grad = constant_factor * euclid_grad;  % 乘以常数constant_factor,防止matlab因为精度问题导致错误
% end

function euclid_grad = euclid_grad_fun(p_R, psi, v_phaseshift, params)
% grad = lambda * theta - beta

[~,~, hc_H, ht_H, ~, Uc, a_bar_BU_H, Ut, a_bar_BT_H, ~, ~] = bulid_H(p_R, psi, v_phaseshift, params);
f_tmp = solve_beamforming_vector(hc_H, ht_H, params);
w = sqrt(params.Pt) * f_tmp;

alpha_B = Uc * w;
alpha_B_tilde = a_bar_BU_H * w;
alpha_E = Ut * w;
alpha_E_tilde = a_bar_BT_H * w;

if isfield(params, 'mu_v2')
    mu = params.mu_v2;
else
    mu = 1;
end

L = length(v_phaseshift);
Phi = alpha_E * alpha_E' - mu * alpha_B * alpha_B';
Phi = (Phi + Phi') / 2;
lambda_max_Phi = max(real(eig(Phi)));
beta = (lambda_max_Phi * eye(L) - Phi) * v_phaseshift + ...
    mu * conj(alpha_B_tilde) * alpha_B - conj(alpha_E_tilde) * alpha_E;

euclid_grad = lambda_max_Phi * v_phaseshift - beta;

if isfield(params, 'grad_clip_norm')
    clip_norm = params.grad_clip_norm;
else
    clip_norm = 1e3;
end
g_norm = norm(euclid_grad, 2);
if g_norm > clip_norm && g_norm > 0
    euclid_grad = euclid_grad * (clip_norm / g_norm);
end
end