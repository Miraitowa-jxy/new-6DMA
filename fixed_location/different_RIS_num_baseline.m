clc;clear all
%% 场景设置
frequency = 3.6e9 ; %载波频率3GHz
c = 3e8; %光速
lambda = c / frequency;
Nt = 32; %BS发射天线数：ULA
Nr = 32; %BS接收天线数：ULA
p_B = [0,0,0].' ; %BS全局坐标
p_U = [280,0,0].' ; %单天线UE全局坐标
p_T = [0, 20, 0].' ; %Target全局坐标
H = 150; %表示RIS的高度一般保持不变
Pt_dB = 20; %发射功率为20dBm
Pt = 1e-3*10^(Pt_dB / 10); 
B = 3e6; % 带宽3MHz
sigma_s2_dB = -174 +10*log10(B); %感知场景噪声方差约为-110dBm
sigma_c2_dB = -174 +10*log10(B); %通信场景噪声方差约为-110dBm
sigma_s2 = 1e-3*10^(sigma_s2_dB / 10);
sigma_c2 = 1e-3*10^(sigma_c2_dB / 10);
gamma_0_dB = 10; %通信场景信噪比阈值10dB
gamma_0 = 10^(gamma_0_dB / 10);
%RIS的行反射单元数目
% Nx_total = [3,4,5,6,7,8,9,10];
Nx_total = [3,4,5,6,7,8,9,10,11,12,13,14];
N_total = Nx_total.^2;
ris_num = length(Nx_total);
constant_factor = 1e20;
% 创建一个结构体来存储参数
params = struct;
params.lambda = lambda;
params.Nt = Nt;
params.Nr = Nr;
params.p_B = p_B;
params.p_U = p_U;
params.p_T = p_T;
params.H = H;
params.Pt = Pt;
params.sigma_s2 = sigma_s2;
params.sigma_c2 = sigma_c2;
params.gamma_0 = gamma_0;
params.constant_factor = constant_factor;
%% 优化问题求解
max_iter = 100; % 最大迭代次数
% 变量个数
num_vars = 3;     
% 定义粒子群大小，变量个数num_vars<= 10时，总数保持在20~50就行;
num_particles = 50; 
p_R_total = [100,100];  %(100,100)作为固定点可以
location_num = size(p_R_total,1);
SNR_s_dB = [];
SNR_c_dB = [];
SNR_c_dB_final = zeros(location_num,ris_num);
SNR_s_dB_final = zeros(location_num,ris_num);
rho_total = [];
rho_total_final = zeros(location_num,ris_num);
normal_vector = cell(ris_num,1);
rho = 0;
aaa_total = [];
aaa = 0;

for j =1:1%location_num %分开跑效果是可以的
    params.p_R = p_R_total(j,:);
    p_R = params.p_R;
    % 设置粒子群的上下界
    x_lb = [0, 0, 0];  % 每个变量的下界
    x_ub = [2*pi, 2*pi, 2*pi]; % 每个变量的上界
    % rng(45)
    % x_opt = x_lb + (x_ub - x_lb) .* rand(1, length(x_lb));
    load x_total.mat
    x_opt = x_total(21,:);
    for k= 1:ris_num
        Nx = Nx_total(k);
        Ny = Nx;
        N = Nx * Ny ;  % RIS的阵元数目
        params.Nx = Nx;
        params.Ny = Ny;
        params.N = N;
        theta = linspace(0,pi,N);
        v_init = exp(1j * theta); %初始化v
        v = v_init';
        psi = x_opt.';

        % 固定 p_R, psi，使用黎曼梯度下降法来优化 v
        [v, f_val_2, info] = rieman_grad_fun_PR(p_R.', psi, x_opt, v, params);
        % 显示当前优化结果
        disp(['第',num2str(j),'个位置,第',num2str(k),'个反射单元数量,信道相关性为：',num2str(rho)]);
 
         % 保存当前值
        psi_final = psi;
        v_final = v; 
 
        [rho,~, hc_H, ht_H, hr, ~, ~, ~, ~, ~, ~] = bulid_H(p_R.', psi_final, v_final, params);
        aaa = sqrt(gamma_0*sigma_c2/Pt/norm(hc_H,2)^2);
        f = solve_beamforming_vector(hc_H, ht_H,params); % 波束成形向量
        % f = hc_H' / norm(hc_H',2);  % 波束成形向量对准感知信道
        SNR_s = Pt * norm(hr*ht_H*f,2)^2/ sigma_s2;
        SNR_c = Pt * norm(hc_H*f,2)^2/ sigma_c2;
        SNR_s_dB_final(j,k) = 10*log10(SNR_s);
        SNR_c_dB_final(j,k) = 10*log10(SNR_c);
        rho_total_final(j,k) =  rho;
             Q = rotation_matrix(psi_final);
        L = Q(:,3);
        normal_vector{k} = L;
    end
end

figure(1)
%左轴
yyaxis left;   
plot(N_total, SNR_s_dB_final(1,:),"o-");hold on;
% plot(N_total, SNR_s_dB_final(2,:),"^-");hold on;
% plot(N_total, SNR_s_dB_final(3,:),"s-");hold on;
ylabel('感知信噪比');

yyaxis right;
plot(N_total, rho_total_final(1,:),"o-");hold on;
% plot(N_total, rho_total_final(2,:),"^-");hold on;
% plot(N_total, rho_total_final(3,:),"s-");hold on;
ylabel('\rho');
legend('location1','location2','location3');
xlabel('RIS反射单元数量');











