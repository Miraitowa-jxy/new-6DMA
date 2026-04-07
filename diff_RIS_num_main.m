clc
clear all
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

p_R_begin = [100,100]; %固定点A=(100,100)
location_num = 3; %第一个为A扩展的小区域
%第二个为A扩展的大区域

% 变量个数
num_vars = 5;     
% 定义粒子群大小，变量个数num_vars<= 10时，总数保持在20~50就行;
num_particles = 200; %对不同区域可以适当修改，区域较小的，50即可，较大的200

limit = zeros(location_num,ris_num);
SNR_c_dB_final = zeros(location_num,ris_num);
SNR_s_dB_final = zeros(location_num,ris_num);
rho_total_final = zeros(location_num,ris_num);
PR_location_end = cell(ris_num,1);
normal_vector = cell(ris_num,1);
rho = 0;

for j=3:3 %一个一个跑效果更好
    if j == 1  % 小区域
        % 设置粒子群的上下界
        x_lb = [50, 50, 0, 0, 0];  % 每个变量的下界
        x_ub = [100, 100, 2*pi, 2*pi, 2*pi]; % 每个变量的上界
        rng(45)
        x_opt = x_lb + (x_ub - x_lb) .* rand(1, length(x_lb));
        x_opt(1:2) = p_R_begin;
        x_opt(3:end) = [3.74, 2, 5.98];
    elseif j==2 % 大区域
        x_lb = [0, 0, 0, 0, 0]; 
        x_ub = [100, 100, 2*pi, 2*pi, 2*pi];
        rng(45)
        x_opt = x_lb + (x_ub - x_lb) .* rand(1, length(x_lb));
        x_opt(1:2) = p_R_begin; 
        x_opt(3:end) = [3.74, 2, 5.98];
    else %固定位置
        x_lb = [100,100,0, 0, 0];  
        x_ub = [100,100,2*pi, 2*pi, 2*pi];
        rng(45)
        x_opt = x_lb + (x_ub - x_lb) .* rand(1, length(x_lb));
        x_opt(3:end) = [3.74, 2, 5.98];
    end  
    for k = 1:ris_num
        Nx = Nx_total(k);
        Ny = Nx;
        N = Nx * Ny ;  % RIS的阵元数目
        params.Nx = Nx;
        params.Ny = Ny;
        params.N = N;
    
        %% 优化问题求解
        max_iter = 100; % 最大迭代次数
        % 初始化v
        theta = linspace(0,pi,N);
        v_init = exp(1j * theta);
        v = v_init';
        vetor_total=[];

        for iter = 1:max_iter
            %定义目标函数
            obj_fun1 = @(x) fixed_v_fun(x, v, params);  % 固定RIS相移向量v之后的目标函数，使用PSO算法求解
        
            % 拉丁超立方（LHS）采样可以提高初始种群的覆盖范围，避免初始点过于集中。
            % rng('shuffle'); 
            initial_swarm = lhsdesign(num_particles, num_vars) .* (x_ub - x_lb) + x_lb;
            
            % 替换第一个粒子为当前最优的x。
            initial_swarm(1, :) = x_opt;
            % 固定 v，使用粒子群优化来优化 p_R 和 psi
            options = optimoptions('particleswarm', ...
            'SwarmSize', num_particles,...             % 变量个数num_vars<= 10时，总数保持在20~50就行;
            'FunctionTolerance', 1e-6     , ...  % 设置目标函数收敛容忍度
            'MaxStallIterations', 50, ...  % 允许优化在目标函数变化很小的情况下继续迭代的最大次数
            'InitialSwarmMatrix', initial_swarm, ...% 'UseParallel', true, ...   % 开启并行计算
            'InertiaRange', [0.1, 0.9], ...% 惯性权重范围
            'SocialAdjustmentWeight', 2, ...       % 适当增加让粒子更倾向于探索全局最优解
            'SelfAdjustmentWeight', 1.6,...%这两个参数选择（1.6，1.6），（1.6，2）
            'MaxIterations', 500); % 最大迭代次数 , ...
            % 'OutputFcn', @pso_output_function); % 绑定迭代输出函数,验证PSO是否是选出最优点         
            [x_opt, f_val_1, exitflag, output] = particleswarm(obj_fun1, num_vars, x_lb, x_ub, options);
            % 更新 p_R, psi
            p_R = x_opt(1:2).';
            psi = x_opt(3:5).';
        
            % 固定 p_R, psi，使用黎曼梯度下降法来优化 v
            [v, f_val_2, info] = rieman_grad_fun(p_R, psi, x_opt, v, params);
        
           % 显示当前优化结果
            disp(['第',num2str(j),'个位置,第',num2str(k),'个反射单元数量,第',num2str(iter),'次迭代, ' 'v更新前函数值: ', ...
                num2str(f_val_1), ', v更新后函数值: ', num2str(f_val_2),',信道相关性为：',num2str(rho)]);
        
             % 保存当前值
            p_R_pre = p_R;
            psi_final = psi;
            v_final = v; 
        
            % 可选：检查是否收敛或达到停止条件
            if iter > max_iter || (abs(f_val_1 - f_val_2) /abs(f_val_1))  < 1e-3
                break;
            end

        end

        [rho,~, hc_H, ht_H, hr, ~, ~, ~, ~, ~, ~] = bulid_H(p_R_pre, psi_final, v_final, params);
        f = solve_beamforming_vector(hc_H, ht_H,params); % 波束成形向量
        limit(j,k) = 10*log10(Pt*norm(hc_H,2)^2/sigma_c2);
        SNR_s = Pt * norm(hr*ht_H*f,2)^2/ sigma_s2;
        SNR_s_dB_final(j,k) = 10*log10(SNR_s);
        SNR_c = Pt * norm(hc_H*f,2)^2/ sigma_c2;
        SNR_c_dB_final(j,k) = 10*log10(SNR_c);
        rho_total_final(j,k) = rho;
        PR_location_end{k} = p_R_pre;
        Q = rotation_matrix(psi_final);
        L = Q(:,3);
        normal_vector{k} = L;
    end
end

figure(1)
%左轴
yyaxis left;   
% plot(N_total, SNR_s_dB_final(1,:),"o-");hold on;
plot(N_total, SNR_s_dB_final(2,:),"^-");hold on;
% plot(N_total, SNR_s_dB_final(3,:),"^-");hold on;
ylabel('SNR of Sensing');
ylim([0 50]);           % 限制横坐标范围
set(gca, 'YTick', 0:5:60); 
yyaxis right;
% plot(N_total, rho_total_final(1,:),"o-");hold on;
plot(N_total, rho_total_final(2,:),"^-");hold on;
% plot(N_total, rho_total_final(3,:),"^-");hold on;
ylabel('\rho');
grid on;
legend('location1,SNR','location1,\rho');
xlabel('Number of RIS reflection units');

  
% figure(1)
% % 左轴
% yyaxis left;   
% plot(N_total, SNR_s_dB_final_AO1,"o-", 'LineWidth', 1.2);hold on;   % 表示小区域的数据
% plot(N_total, SNR_s_dB_final_AO2,"^-", 'LineWidth', 1.2);hold on; % 表示大区域的数据
% plot(N_total, SNR_s_dB_final_wA,"*:", 'LineWidth', 1.2);hold on;  % 表示固定pos的数据,但改变RIS姿态
% plot(N_total, SNR_s_dB_final_baseline,"s:", 'LineWidth', 1.2);hold on;%对比方法，只改变相移向量v，不改变RIS的位置和姿态
% ylabel('SNR(dB) of Sensing','FontSize', 14);
% ylim([10 55]);           % 限制横坐标范围
% yticks(10:5:55);        % 设置横坐标刻度 
% xticks(0:20:200);
% yyaxis right;
% plot(N_total, rho_total_final_AO1,"o-", 'LineWidth', 1.2);hold on;
% plot(N_total, rho_total_final_AO2,"^-", 'LineWidth', 1.2);hold on;
% plot(N_total, rho_total_final_wA,"*:", 'LineWidth', 1.2);hold on;
% plot(N_total, rho_total_final_baseline,"s:", 'LineWidth', 1.2);hold on;
% ylabel('\rho','FontSize', 14);
% grid on;
% h_legend = legend('6D + PBF (R1),SNR','6D + PBF (R2),SNR','3D orientation + PBF,SNR','PBF Only,SNR',...
%      '6D + PBF (R1),\rho','6D + PBF (R2),\rho','3D orientation + PBF,\rho','PBF Only,\rho','FontSize', 11);
% set(h_legend, 'LineWidth', 1.5); % 图例加粗
% set(gca, 'LineWidth', 1.5,'FontSize', 12);  % 加粗坐标轴
% xlabel('Number of reflecting elements','FontSize', 14);


