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

function euclid_grad = euclid_grad_fun(p_R, psi, theta_phase, params)
%EUCLID_GRAD_FUN Euclidean gradient of the V2 (10) phase surrogate.
% The manifold variable here is the paper variable theta = conj(diag(Phi)).
% The channel builder still expects the physical reflection vector phi=diag(Phi),
% so this function converts theta back to phi=conj(theta) only for channel terms.

channel_params = params;
channel_params.skip_objective_eval = true;
phi_for_channel = conj(theta_phase);
[~,~, hc_H, ht_H, ~, Uc, a_bar_BU_H, Ut, a_bar_BT_H, ~, ~] = ...
    bulid_H(p_R, psi, phi_for_channel, channel_params);

[lambda_max_Phi, beta, ~] = v2_surrogate_terms(theta_phase, hc_H, ht_H, Uc, a_bar_BU_H, Ut, a_bar_BT_H, params);

% Teacher-specified phase gradient. It is one half of the exact Euclidean
% gradient of Eq. (10), i.e., 0.5 * (2*lambda*theta - 2*beta), so it keeps
% the same stationary point and descent direction while following the formula
% requested by the advisor.
euclid_grad = lambda_max_Phi * theta_phase - beta;
if isfield(params, 'phase_objective_scale')
    euclid_grad = params.phase_objective_scale * euclid_grad;
elseif isfield(params, 'constant_factor')
    euclid_grad = params.constant_factor * euclid_grad;
end

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
