function euclid_grad = euclid_grad_fun(p_R, psi, v_phaseshift, params)
[~,~, hc_H, ht_H, hr, Uc, a_bar_BU_H, Ut, a_bar_BT_H, Ur, a_bar_TB] = bulid_H(p_R, psi, v_phaseshift, params);
%将目标函数分为三个部分
f1 = hc_H * ht_H';
f2 = -hr' * hr;
f3 = ht_H * hc_H';
%对各部分函数求导
delta_f1 = conj(Ut) * a_bar_BU_H.' + conj(Ut) * Uc.' * v_phaseshift;  
delta_f2 = -Ur' * a_bar_TB - Ur' * Ur * v_phaseshift;
delta_f3 = conj(Uc) * a_bar_BT_H.' + conj(Uc) * Ut.' * v_phaseshift;
euclid_grad = f2*f3*delta_f1 + f3*f1*delta_f2 + f1*f2*delta_f3;
constant_factor = params.constant_factor;
euclid_grad = constant_factor * euclid_grad;  % 乘以常数constant_factor,防止matlab因为精度问题导致错误
end