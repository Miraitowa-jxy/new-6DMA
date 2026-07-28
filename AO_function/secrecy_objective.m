function [result, secrecy_rate, snr_B, snr_E] = secrecy_objective(theta, hc_H, ht_H, Uc, a_bar_BU_H, Ut, a_bar_BT_H, params)

[lambda_max_Phi, beta, f_tmp] = v2_surrogate_terms(theta, hc_H, ht_H, Uc, a_bar_BU_H, Ut, a_bar_BT_H, params);
% V2 paper Eq. (10) phase surrogate objective (constant term c omitted):
%   minimize lambda_max(Phi) * ||theta||^2 - 2 * real(theta' * beta),
%   subject to |theta_i| = 1.
% This is not the original secrecy-rate objective; it is the MM/bounding
% subproblem used to update the IRS phase vector with fixed beamformer.
result = lambda_max_Phi * real(theta' * theta) - 2 * real(theta' * beta);

snr_B = params.Pt * norm(hc_H * f_tmp, 2)^2 / max(params.sigma_c2, 1e-20);
snr_E = params.Pt * norm(ht_H * f_tmp, 2)^2 / max(params.sigma_s2, 1e-20);
secrecy_rate = log2((1 + snr_B) / max(1 + snr_E, 1e-20));
end
