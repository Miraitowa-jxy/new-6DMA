%% ISAC信道建模
function [rho,result, hc_H, ht_H, hr, Uc, a_bar_BU_H, Ut, a_bar_BT_H, Ur, a_bar_TB] = bulid_H(p_R_init, psi, v_phaseshift, params)
%固定参数
lambda = params.lambda;
Nt = params.Nt;
Nr = params.Nr;
Nx = params.Nx;
Ny = params.Ny;
p_B = params.p_B;
p_U = params.p_U;
p_T = params.p_T;
H = params.H; %表示RIS的高度一般保持不变
%通信场景:hc_H,感知场景:ht_H,hr
dB = lambda / 2; %BS天线间距
dR = lambda / 2; %RIS天线间距
p_R = [p_R_init;H]; %RIS左下角元素全局坐标
PHI = diag(v_phaseshift);
Q = rotation_matrix(psi);
p_BL = Q.' * (p_B - p_R); %BS局部坐标
p_UL = Q.' * (p_U - p_R); %单天线UE局部坐标
p_TL = Q.' * (p_T - p_R); %Target局部坐标
p_RL = Q.' * (p_R - p_R); %Target局部坐标

% ULA阵列响应向量
phi_BR = acos((p_R(3) - p_B(3))/norm(p_R - p_B, 2));
a_BR = array_response_ULA(phi_BR, Nt, dB, lambda);

phi_RB = acos((p_B(3) - p_R(3))/norm(p_B - p_R, 2));
a_RB = array_response_ULA(phi_RB, Nr, dB, lambda);

phi_BU = acos((p_U(3) - p_B(3))/norm(p_U - p_B, 2));
a_BU = array_response_ULA(phi_BU, Nt, dB, lambda);

phi_BT = acos((p_T(3) - p_B(3))/norm(p_T - p_B, 2));
a_BT = array_response_ULA(phi_BT, Nt, dB, lambda);

phi_TB = acos((p_B(3) - p_T(3))/norm(p_B - p_T, 2));
a_TB = array_response_ULA(phi_TB, Nr, dB, lambda);

% UPA阵列响应向量
q_RB = p_RL - p_BL;
theta_e_RBA = acos(q_RB(3)/norm(q_RB, 2));
theta_a_RBA = atan2(q_RB(2),q_RB(1));
a_A = array_response_UPA(theta_e_RBA, theta_a_RBA, Nx, Ny, dR, lambda);

q_BR = p_BL - p_RL;
theta_e_RBD = acos(q_BR(3)/norm(q_BR, 2));
theta_a_RBD = atan2(q_BR(2),q_BR(1));
a_B = array_response_UPA(theta_e_RBD, theta_a_RBD, Nx, Ny, dR, lambda);

q_UR = p_UL - p_RL;
theta_e_RUD = acos(q_UR(3)/norm(q_UR, 2));
theta_a_RUD = atan2(q_UR(2),q_UR(1));
a_DU = array_response_UPA(theta_e_RUD, theta_a_RUD, Nx, Ny, dR, lambda);

q_TR = p_TL - p_RL;
theta_e_RTD = acos(q_TR(3)/norm(q_TR, 2));
theta_a_RTD = atan2(q_TR(2),q_TR(1));
a_DT = array_response_UPA(theta_e_RTD, theta_a_RTD, Nx, Ny, dR, lambda);

q_RT = p_RL - p_TL;
theta_e_RTA = acos(q_RT(3)/norm(q_RT, 2));
theta_a_RTA = atan2(q_RT(2),q_RT(1));
a_AT = array_response_UPA(theta_e_RTA, theta_a_RTA, Nx, Ny, dR, lambda);

% 复信道系数
% Gt_dBi = 10*log10(Nt);
% Gr_dBi = 10*log10(Nr);
Gt_dBi = 0;
Gr_dBi = 0;
Gt = 10^(Gt_dBi / 10);
Gr = 10^(Gr_dBi / 10);
Gr_UE = 1;

% % 初始化相位
% phase = [0.37,0.42,0.89,0.14,0.91,0.56];
% factor1 = 3;
% factor2 = 2.2;
% alpha_cLOS = sqrt((lambda/4/pi)^2 * Gt * Gr_UE* (1/norm(p_U -p_B,2))^factor1) * exp(1j*2 * pi * phase(1));
% alpha_cNLOS = sqrt((lambda/4/pi)^2 * Gt * Gr_UE * (1/(norm(p_U -p_R,2)+norm(p_R -p_B,2)))^factor1) * exp(1j*2 * pi * phase(2));
% alpha_tLOS = sqrt((lambda/4/pi)^2 * Gt * (1/norm(p_T -p_B,2))^factor2) * exp(1j*2 * pi * phase(3));
% alpha_tNLOS = sqrt((lambda/4/pi)^2 * Gt *(1/(norm(p_T -p_R,2)+norm(p_R -p_B,2)))^factor2) * exp(1j*2 * pi * phase(4));
% alpha_rLOS = sqrt((lambda/4/pi)^2 * Gr * (1/norm(p_T -p_B,2))^factor2) * exp(1j*2 * pi * phase(5));
% alpha_rNLOS = sqrt((lambda/4/pi)^2 * Gr * (1/(norm(p_T -p_R,2)+norm(p_R -p_B,2)))^factor2) * exp(1j*2 * pi * phase(6));

% 初始化相位
factor1 = 3;
factor2 = 2.2;
alpha_cLOS = sqrt((lambda/4/pi)^2 * Gt * Gr_UE* (1/norm(p_U -p_B,2))^factor1) * exp(1j*2 * pi * norm(p_U -p_B,2) / lambda);
alpha_cNLOS = sqrt((lambda/4/pi)^2 * Gt * Gr_UE * (1/(norm(p_U -p_R,2)+norm(p_R -p_B,2)))^factor1) * exp(1j*2 * pi * (norm(p_U -p_R,2)+norm(p_R -p_B,2))/ lambda);
alpha_tLOS = sqrt((lambda/4/pi)^2 * Gt * (1/norm(p_T -p_B,2))^factor2) * exp(1j*2 * pi * norm(p_T -p_B,2)/ lambda);
alpha_tNLOS = sqrt((lambda/4/pi)^2 * Gt *(1/(norm(p_T -p_R,2)+norm(p_R -p_B,2)))^factor2) * exp(1j*2 * pi * (norm(p_T -p_R,2)+norm(p_R -p_B,2))/ lambda);
alpha_rLOS = sqrt((lambda/4/pi)^2 * Gr * (1/norm(p_T -p_B,2))^factor2) * exp(1j*2 * pi * norm(p_T -p_B,2)/ lambda);
alpha_rNLOS = sqrt((lambda/4/pi)^2 * Gr * (1/(norm(p_T -p_R,2)+norm(p_R -p_B,2)))^factor2) * exp(1j*2 * pi * (norm(p_T -p_R,2)+norm(p_R -p_B,2))/ lambda);

%% 有效孔径增益
L = Q(:,3);
norm_L = norm(L,2);
cos_phi1 = ((p_R - p_B).' * L) / (norm(p_R - p_B,2) * norm_L);
cos_phi2 = ((p_R - p_U).' * L) / (norm(p_R - p_U,2) * norm_L);
cos_phi3 = ((p_R - p_T).' * L) / (norm(p_R - p_T,2) * norm_L);
F_RU = cos_phi1 * cos_phi2;
F_RT = cos_phi1 * cos_phi3;

%% 信道建模
%通信场景
hc_H = sqrt(F_RU)*(alpha_cLOS * a_BU' + alpha_cNLOS * (a_DU' * PHI * a_A) * a_BR'); % 这是个标量
%感知场景
ht_H = sqrt(F_RT)*(alpha_tLOS * a_BT' + alpha_tNLOS  *  a_DT' * PHI * a_A * a_BR');
hr = sqrt(F_RT)*(alpha_rLOS * a_TB + alpha_rNLOS * a_RB * a_B' * PHI * a_AT);
result = -norm(hr*ht_H*hc_H',2)^2;

rho = abs(hc_H*ht_H')/norm(hc_H,2)/norm(ht_H,2);
%% 简化的表示
Uc = sqrt(F_RU)*alpha_cNLOS * diag(a_DU') * a_A * a_BR';
a_bar_BU_H = sqrt(F_RU)*alpha_cLOS * a_BU';
Ut = sqrt(F_RT)*alpha_tNLOS  *  diag(a_DT') * a_A * a_BR';
a_bar_BT_H = sqrt(F_RT)*alpha_tLOS * a_BT';
Ur = sqrt(F_RT)*alpha_rNLOS * a_RB * a_B' * diag(a_AT);
a_bar_TB = sqrt(F_RT)*alpha_rLOS * a_TB;
end

function a_R = array_response_ULA(phi, N, d, lambda)
    % phi: 俯仰角
    % N: 天线数量
    % d: 天线间距
    % lambda: 波长
    n = 0:N-1;  
    % 生成向量
    a_R = exp(1j * 2 * pi / lambda * d * cos(phi) * n).';  % .'转置
end

function a_R = array_response_UPA(theta_elevation, theta_azimuth, Nx, Ny, d, lambda)
    % phi_elevation: 俯仰角
    % phi_azimuth: 方位角
    % Nx: 水平天线数量
    % Ny: 垂直天线数量
    % d: 天线间距
    % lambda: 波长
    n = 0:Nx-1;
    m = 0:Ny-1;
    % 生成阵列相应向量
    a_x = exp(1j * 2 * pi / lambda * d * sin(theta_elevation) * cos(theta_azimuth) * n).';  % 列向量
    a_y = exp(1j * 2 * pi / lambda * d * sin(theta_elevation) * sin(theta_azimuth) * m).';
    a_R = kron(a_x,a_y);
end