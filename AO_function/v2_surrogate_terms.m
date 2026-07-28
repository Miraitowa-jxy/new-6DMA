function [lambda_max_Phi, beta, f_tmp] = v2_surrogate_terms(theta, hc_H, ht_H, Uc, a_bar_BU_H, Ut, a_bar_BT_H, params)
%V2_SURROGATE_TERMS Shared lambda/beta construction for the V2 (10) phase block.
% Keep theta_prev and the beamformer fixed during one AO phase subproblem so
% secrecy_objective and euclid_grad_fun use identical surrogate parameters.

if isfield(params, 'fixed_beamformer') && ~isempty(params.fixed_beamformer)
    f_tmp = params.fixed_beamformer;
else
    f_tmp = solve_beamforming_vector(hc_H, ht_H, params);
end
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
if isfield(params, 'theta_prev') && ~isempty(params.theta_prev)
    theta_prev = params.theta_prev;
else
    theta_prev = theta;
end

L = length(theta);
Phi = alpha_E * alpha_E' - mu * alpha_B * alpha_B';
Phi = (Phi + Phi') / 2;
lambda_max_Phi = max(real(eig(Phi)));
beta = (lambda_max_Phi * eye(L) - Phi) * theta_prev + ...
    mu * conj(alpha_B_tilde) * alpha_B - conj(alpha_E_tilde) * alpha_E;
end