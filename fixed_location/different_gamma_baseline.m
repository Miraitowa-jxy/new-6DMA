clc;clear all
%% 体现tradeoff特性，基线方案
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

gamma_0_total =25:0.1:31.7; %通信场景信噪比阈值
gamma_num = length(gamma_0_total);

%RIS的行反射单元数目
Nx = 4;
Ny = Nx;
N = Nx * Ny ;  % RIS的阵元数目
constant_factor = 1e20;
% 创建一个结构体来存储参数
params = struct;
params.lambda = lambda;
params.Nt = Nt;
params.Nr = Nr;
params.Nx = Nx;
params.Ny = Ny;
params.N = N;
params.p_B = p_B;
params.p_U = p_U;
params.p_T = p_T;
params.H = H;
params.Pt = Pt;
params.sigma_s2 = sigma_s2;
params.sigma_c2 = sigma_c2;
params.constant_factor = constant_factor;
%% 优化问题求解
max_iter = 100; % 最大迭代次数 
 
p_R = [100,100]; % 
params.p_R = p_R;
location_num = 1000;
SNR_c_total_final = zeros(location_num,gamma_num);
SNR_s_total_final = zeros(location_num,gamma_num);
compare_rho = 0;
rho_total = zeros(location_num,1);
rho_total_final = zeros(location_num,gamma_num);
rho = 0;
aaa_total = [];
% load v.mat

x_total = zeros(location_num,3);   
% 设置粒子群的上下界
x_lb = [0, 0, 0];  % 每个变量的下界
x_ub = [2*pi, 2*pi, 2*pi]; % 每个变量的上界
% for num =1:location_num
%     x_num = x_lb + (x_ub - x_lb) .* rand(1, length(x_lb));
%     x_total(num,:) = x_num;
% end
load x_total.mat
limit = zeros(location_num,1);
rho_total = [];
compare_rho_total = [];
for j =1:1 %location_num
    x_opt = x_total(21,:);%21
    for k= 1:gamma_num   
        gamma_0 = 10^(gamma_0_total(k) / 10);%通信场景信噪比阈值
        params.gamma_0 = gamma_0;
        psi = x_opt.';
        %初始化v
        theta = linspace(0,pi,N);
        v_init = exp(1j * theta); 
        v = v_init';
        % 固定 p_R, psi，使用黎曼梯度下降法来优化 v
        [v, f_val_2, info] = rieman_grad_fun_PR(p_R.', psi, x_opt, v, params);
 
         % 保存当前值
        psi_final = psi;
        v_final = v; 

        [rho,~, hc_H, ht_H, hr, ~, ~, ~, ~, ~, ~] = bulid_H(p_R.', psi_final, v_final, params);
        limit(j,:) = 10*log10(Pt*norm(hc_H,2)^2/sigma_c2);
        % %求出在给定hc的情况下，通信阈值的上界，必须严格小于等于这个值。
        aaa = sqrt(gamma_0*sigma_c2/Pt/norm(hc_H,2)^2);
        rho_total(j,1) = rho;
        aaa_total = [aaa_total,aaa];
        f = solve_beamforming_vector(hc_H, ht_H,params); % 波束成形向量
        % f = hc_H' / norm(hc_H',2); 
         disp(['第',num2str(j),'个姿态,第',num2str(k),'个通信阈值,信道相关性为：',num2str(rho)]);
        SNR_s = Pt * norm(hr*ht_H*f,2)^2/ sigma_s2;
        SNR_c = Pt * norm(hc_H*f,2)^2/ sigma_c2;

        SNR_s_total_final(j,k) = 10*log10(SNR_s);
        SNR_c_total_final(j,k) = 10*log10(SNR_c);
     
    end
end

x1 = SNR_c_total_final(1,:);
y1 = SNR_s_total_final(1,:);%固定位置1，但优化姿态
figure(1) 
plot(x1, y1,"-", 'LineWidth', 2);hold on;
grid on;
xlabel('SNR(dB) of Communication');
ylabel('SNR(dB) of Sensing');
legend('fixed-pos1-wA(\rho = 0.575)');

% figure(1) 
% plot(SNR_c_total_final_baseline, SNR_s_total_final_baseline,"-", 'LineWidth', 2);hold on;
% grid on;
% xlabel('SNR(dB) of Communication');
% ylabel('SNR(dB) of Sensing');
% legend('fixed-pos1-wA(\rho = 0.67)');