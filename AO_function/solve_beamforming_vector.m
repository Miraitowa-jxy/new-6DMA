% function f = solve_beamforming_vector(hc_H, ht_H,params)
% hc = hc_H';
% ht = ht_H';
% Pt = params.Pt;
% % sigma_s2 = params.sigma_s2;
% sigma_c2 = params.sigma_c2;
% gamma_0 = params.gamma_0;
% 
% %利用文章《RIS-assisted Integrated Sensing and Communications: A Subspace
% %Rotation Approach》的定理一求解
% condition_1 = Pt*abs(hc_H*ht)^2;
% condition_2 = gamma_0*sigma_c2*norm(ht,2)^2;
% if condition_1 >= condition_2
%     f = ht / norm(ht,2);
% else
%     u1 = hc / norm(hc,2); 
%     u2 = (ht - (u1' * ht) * u1) / norm(ht - (u1' * ht) * u1,2);  
%     x11 =sqrt(gamma_0*sigma_c2/norm(hc,2)^2) * (u1' * ht) / abs(u1' * ht);
%     x1 = x11 /sqrt(Pt);
%     x22 =sqrt(Pt-gamma_0*sigma_c2/norm(hc,2)^2) * (u2' * ht) / abs(u2' * ht);
%     x2 = x22 /sqrt(Pt);
%     f = x1*u1 +x2*u2;
% end 
% end

function f = solve_beamforming_vector(hB_H, hE_H, params)

hB = hB_H';
hE = hE_H';
Nt = length(hB);

if norm(hB, 2) <= 1e-12
    f = zeros(size(hB));
    return;
end

A = eye(Nt) + (params.Pt / max(params.sigma_c2, 1e-20)) * (hB * hB');
B = eye(Nt) + (params.Pt / max(params.sigma_s2, 1e-20)) * (hE * hE');
A = (A + A') / 2;
B = (B + B') / 2;

[V, D] = eig(A, B, 'vector');
[~, idx] = max(real(D));
f = V(:, idx);
f = f / max(norm(f, 2), 1e-12);

end
