function [cost, secrecy_rate, snr_B, snr_E, f_tmp] = secrecy_rate_objective(hB_H, hE_H, params)
%SECRECY_RATE_OBJECTIVE True V2 Eq. (3)-style objective for fixed g/Phi/f.
% In the 6DMA model, hB_H = h_B(g,Phi)^H and hE_H = h_E(g,Phi)^H.
% We minimize the negative secrecy rate so it can be used by PSO/Manopt-style
% minimizers while the physical objective remains secrecy-rate maximization.

if isfield(params, 'fixed_beamformer') && ~isempty(params.fixed_beamformer)
    f_tmp = params.fixed_beamformer;
else
    f_tmp = solve_beamforming_vector(hB_H, hE_H, params);
end

snr_B = params.Pt * abs(hB_H * f_tmp)^2 / max(params.sigma_c2, 1e-20);
snr_E = params.Pt * abs(hE_H * f_tmp)^2 / max(params.sigma_s2, 1e-20);
% Keep the signed secrecy rate during optimization. Clipping negative
% values to zero makes many candidate g/Phi/f points indistinguishable and
% can falsely look like instant AO convergence.
secrecy_rate = log2((1 + snr_B) / max(1 + snr_E, 1e-20));
cost = -secrecy_rate;
end