function [result, secrecy_rate, snr_B, snr_E] = secrecy_objective(theta, hc_H, ht_H, Uc, a_bar_BU_H, Ut, a_bar_BT_H, params)

if isfield(params, 'fixed_beamformer') && ~isempty(params.fixed_beamformer)
    f_tmp = params.fixed_beamformer;
else
    f_tmp = solve_beamforming_vector(hc_H, ht_H, params);
end
w = sqrt(params.Pt) * f_tmp;
L = length(theta);

alpha_B = Uc * w;
alpha_B_tilde = a_bar_BU_H * w;
alpha_E = Ut * w;
alpha_E_tilde = a_bar_BT_H * w;

if isfield(params, 'mu_v2')
    mu = params.mu_v2;
else
    mu = 1;
end
if isfield(params, 'theta_prev')
    theta_prev = params.theta_prev;
else
    theta_prev = theta;
end

Phi = alpha_E * alpha_E' - mu * alpha_B * alpha_B';
Phi = (Phi + Phi') / 2;
lambda_max_Phi = max(real(eig(Phi)));
beta = (lambda_max_Phi * eye(L) - Phi) * theta_prev + ...
    mu * conj(alpha_B_tilde) * alpha_B - conj(alpha_E_tilde) * alpha_E;
result = lambda_max_Phi * L - 2 * real(theta' * beta);

snr_B = params.Pt * norm(hc_H * f_tmp, 2)^2 / max(params.sigma_c2, 1e-20);
snr_E = params.Pt * norm(ht_H * f_tmp, 2)^2 / max(params.sigma_s2, 1e-20);
secrecy_rate = max(log2((1 + snr_B) / max(1 + snr_E, 1e-20)), 0);
end