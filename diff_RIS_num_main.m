clc
clear all
root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'AO_function'), '-begin');
addpath(fullfile(root_dir, 'fixed_location'), '-end');
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
constant_factor = 1; % 数值稳定：避免目标/梯度被过度放大导致优化停滞

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
params.secrecy_rate_min = 0; % v2公式10约束阈值（默认无硬门限）
params.secrecy_leakage_weight = 0.5; % 波束对"E端”泄漏抑制权重
params.objective_type = 'v2_eq10_surrogate';
params.mu_v2 = 1;
params.grad_eps_phase = 1e-4;
params.manopt_verbosity = 0;
params.manopt_maxiter = 120;
params.manopt_tolgradnorm = 1e-5;
params.manopt_minstepsize = 1e-10;
params.ao_rel_gap_tol = 1e-3;        % AO块间目标差收敛阈值
params.ao_rel_progress_tol = 5e-4;   % AO整轮目标相对变化收敛阈值

p_R_begin = [100,100]; %固定点A=(100,100)
location_num = 3; %第一个为A扩展的小区域
%第二个为A扩展的大区域

% 变量个数
num_vars = 5;     
% 定义粒子群大小，变量个数num_vars<= 10时，总数保持在20~50就行;
num_particles = 50; %对不同区域可以适当修改，区域较小的，50即可，较大的200

limit = zeros(location_num,ris_num);
SNR_c_dB_final = zeros(location_num,ris_num);
SNR_s_dB_final = zeros(location_num,ris_num);
rho_total_final = zeros(location_num,ris_num);
PR_location_end = cell(ris_num,1);
normal_vector = cell(ris_num,1);
convergence_history = cell(location_num, ris_num);
g_obj_history = cell(location_num, ris_num);     % 6D参数块（g）更新后的目标值
phi_obj_history = cell(location_num, ris_num);   % 相位块（Phi/θ）更新后的目标值
f_obj_history = cell(location_num, ris_num);     % 波束块（f）对应目标值
loss_history = cell(location_num, ris_num);      % 总体loss
total_obj_history = cell(location_num, ris_num); % 每轮真实总目标（当前x,v）
g_value_history = cell(location_num, ris_num);     % 每轮g块变量值：[p_Rx,p_Ry,psi1,psi2,psi3]
phi_value_history = cell(location_num, ris_num);   % 每轮Phi块变量值：diag(Phi)=v
f_value_history = cell(location_num, ris_num);     % 每轮f块变量值：波束向量f
g_delta_history = cell(location_num, ris_num);     % 每轮g块变量更新范数
phi_delta_history = cell(location_num, ris_num);   % 每轮Phi块变量更新范数
f_delta_history = cell(location_num, ris_num);     % 每轮f块变量更新范数
rho = 0;

% 是否固定一个N做收敛剖析
run_single_N = true;
target_k_track = 6; % 对应 Nx_total(6)=8, N=64

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
    if run_single_N
        target_k_track = min(max(target_k_track,1),ris_num);
        k_range = target_k_track;
    else
        k_range = 1:ris_num;
    end
    for k = k_range
        Nx = Nx_total(k);
        Ny = Nx;
        N = Nx * Ny ;  % RIS的阵元数目
        params.Nx = Nx;
        params.Ny = Ny;
        params.N = N;
    
        %% 优化问题求解
        max_iter = 1e20; % 最大迭代次数上限，实际停止由AO收敛条件决定
        % 初始化v
        theta = linspace(0,pi,N);
        v_init = exp(1j * theta);
        v = v_init';
        vetor_total=[];

        init_p_R = x_opt(1:2).';
        init_psi = x_opt(3:5).';
        [~,~, hc_init, ht_init, ~, ~, ~, ~, ~, ~, ~] = bulid_H(init_p_R, init_psi, v, params);
        f_current = solve_beamforming_vector(hc_init, ht_init, params);

        f_val_2_pre = inf;
        f_curve = [];
        g_curve = [];
        phi_curve = [];
        f_curve_beam = [];
        loss_curve = [];
        g_value_curve = [];
        phi_value_curve = [];
        f_value_curve = [];
        g_delta_curve = [];
        phi_delta_curve = [];
        f_delta_curve = [];
        x_prev_for_delta = x_opt;
        v_prev_for_delta = v;
        f_prev_for_delta = f_current;
        for iter = 1:max_iter
            params.theta_prev = v;
            % g块和Phi块优化时固定上一轮波束，f块在本轮最后单独更新
            params.fixed_beamformer = f_current;
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
            'Display', 'off', ...
            'FunctionTolerance', 1e-4     , ...  % 设置目标函数收敛容忍度
            'MaxStallIterations', 20, ...  % 允许优化在目标函数变化很小的情况下继续迭代的最大次数
            'InitialSwarmMatrix', initial_swarm, ...% 'UseParallel', true, ...   % 开启并行计算
            'InertiaRange', [0.1, 0.9], ...% 惯性权重范围
            'SocialAdjustmentWeight', 2, ...       % 适当增加让粒子更倾向于探索全局最优解
            'SelfAdjustmentWeight', 1.6,...%这两个参数选择（1.6，1.6），（1.6，2）
            'MaxIterations', 200); % 最大迭代次数 , ...
            % 'OutputFcn', @pso_output_function); % 绑定迭代输出函数,验证PSO是否是选出最优点         
            [x_opt, f_val_1, exitflag, output] = particleswarm(obj_fun1, num_vars, x_lb, x_ub, options);
            g_value_curve = [g_value_curve; x_opt];
            g_delta_curve = [g_delta_curve; norm(x_opt - x_prev_for_delta, 2)];
            % 更新 p_R, psi
            p_R = x_opt(1:2).';
            psi = x_opt(3:5).';
        
            % 固定 p_R, psi，使用黎曼梯度下降法来优化 v
            [v, f_val_2, info] = rieman_grad_fun(p_R, psi, x_opt, v, params);
            g_curve = [g_curve; f_val_1];
            phi_curve = [phi_curve; f_val_2];
            phi_value_curve = [phi_value_curve; v(:).'];
            phi_delta_curve = [phi_delta_curve; norm(v - v_prev_for_delta, 2)];

            % f块：固定当前g和Phi后，重新构造/更新波束，并记录更新后的真实目标
            [rho_iter,~, hc_tmp, ht_tmp, ~, Uc_tmp, a_BU_tmp, Ut_tmp, a_BT_tmp, ~, ~] = bulid_H(p_R, psi, v, params);
            f_tmp = solve_beamforming_vector(hc_tmp, ht_tmp, params);
            params.fixed_beamformer = f_tmp;
            [f_obj,~,~,~] = secrecy_objective(v, hc_tmp, ht_tmp, Uc_tmp, a_BU_tmp, Ut_tmp, a_BT_tmp, params);
            f_curve_beam = [f_curve_beam; f_obj];
            loss_curve = [loss_curve; f_obj];
            f_curve = [f_curve; f_obj];
            f_value_curve = [f_value_curve; f_tmp(:).'];
            f_delta_curve = [f_delta_curve; norm(f_tmp - f_prev_for_delta, 2)];

            % 下一轮g/Phi块固定本轮刚更新的波束
            f_current = f_tmp;
            x_prev_for_delta = x_opt;
            v_prev_for_delta = v;
            f_prev_for_delta = f_tmp;
            params = rmfield(params, 'fixed_beamformer');
        
           % 显示当前优化结果
            disp(['第',num2str(j),'个位置,第',num2str(k),'个反射单元数量,第',num2str(iter),'次迭代, ' 'v更新前函数值: ', ...
                num2str(f_val_1), ', Phi更新后函数值: ', num2str(f_val_2), ', f更新后函数值: ', num2str(f_obj), ',信道相关性为：',num2str(rho_iter)]);
            rho = rho_iter;
        
             % 保存当前值
            p_R_pre = p_R;
            psi_final = psi;
            v_final = v; 
        
            % 可选：检查是否收敛或达到停止条件
            rel_gap = abs(f_val_1 - f_obj) / max(abs(f_val_1), 1e-12);
            if isfinite(f_val_2_pre)
                rel_progress = abs(f_obj - f_val_2_pre) / max(abs(f_val_2_pre), 1e-12);
            else
                rel_progress = inf;
            end
            if rel_gap < params.ao_rel_gap_tol || rel_progress < params.ao_rel_progress_tol
                break;
            end
            f_val_2_pre = f_obj;

        end
        if isempty(f_curve)
            f_curve = f_val_2_pre;
        end
        convergence_history{j, k} = f_curve;
        g_obj_history{j, k} = g_curve;
        phi_obj_history{j, k} = phi_curve;
        f_obj_history{j, k} = f_curve_beam;
        loss_history{j, k} = loss_curve;
        total_obj_history{j, k} = loss_curve;
        g_value_history{j, k} = g_value_curve;
        phi_value_history{j, k} = phi_value_curve;
        f_value_history{j, k} = f_value_curve;
        g_delta_history{j, k} = g_delta_curve;
        phi_delta_history{j, k} = phi_delta_curve;
        f_delta_history{j, k} = f_delta_curve;

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

% 收敛折线图（示例：当前运行的 j=3，k=4）
figure(2);
nonempty_mask = ~cellfun(@isempty, convergence_history);
[j_idx, k_idx] = find(nonempty_mask);
if ~isempty(j_idx)
    % 自动选取最后一个有数据的(j,k)，避免硬编码索引没有数据
    target_j = j_idx(end);
    target_k = k_idx(end);
    plot(1:length(convergence_history{target_j, target_k}), convergence_history{target_j, target_k}, 'o-', 'LineWidth', 1.2);
    xlabel('AO Iteration');
    ylabel('Objective Value (after v update)');
    title(['Convergence Curve (j=', num2str(target_j), ', k=', num2str(target_k), ')']);
    grid on;
    % 图2-2：相邻迭代目标变化（观察是否进入平台期）
    if length(convergence_history{target_j, target_k}) >= 2
        figure(3);
        delta_obj = abs(diff(convergence_history{target_j, target_k}));
        semilogy(2:length(convergence_history{target_j, target_k}), delta_obj, 's-', 'LineWidth', 1.1);
        xlabel('AO Iteration');
        ylabel('|f_t - f_{t-1}| (log scale)');
        title(['Convergence gap (j=', num2str(target_j), ', k=', num2str(target_k), ')']);
        grid on;
    end

    % 图2-3：同一位置下，不同RIS规模的收敛曲线对比（最多4条）。
    % run_single_N=true时只会有一个k有数据，因此只画实际存在的数据，避免空白图。
    figure(4);
    hold on;
    available_k = find(~cellfun(@isempty, convergence_history(target_j, :)));
    k_show = available_k(1:min(4, numel(available_k)));
    legend_str = {};
    for kk = k_show
        plot(convergence_history{target_j, kk}, 'LineWidth', 1.1);
        legend_str{end+1} = ['N=', num2str((Nx_total(kk))^2)]; %#ok<SAGROW>
    end
    xlabel('AO Iteration');
    ylabel('Objective Value');
    title(['Convergence comparison across available RIS sizes (j=', num2str(target_j), ')']);
    grid on;
    if ~isempty(legend_str)
        legend(legend_str, 'Location', 'best');
    else
        text(0.1,0.5,'No RIS-size comparison data available.'); axis off;
    end
    hold off;
else
    text(0.1, 0.5, 'No convergence data generated in this run.');
    axis off;
end

% 固定N时，画三个块变量自身的变化；AO总目标只单独画一张图
if run_single_N && ~isempty(g_value_history{j, k_range(1)})
    jj = j;
    kk = k_range(1);
    g_track = g_value_history{jj,kk};
    phi_track = phi_value_history{jj,kk};
    f_track = f_value_history{jj,kk};
    loss_track = total_obj_history{jj,kk};

    figure(5);
    subplot(2,1,1);
    plot(g_track(:,1:2), 'LineWidth', 1.2);
    xlabel('AO Iteration'); ylabel('RIS position / m');
    title(['g block position trajectory (N=', num2str(Nx_total(kk)^2), ')']);
    legend('p_R_x','p_R_y','Location','best');
    grid on;

    subplot(2,1,2);
    plot(g_track(:,3:5), 'LineWidth', 1.2);
    xlabel('AO Iteration'); ylabel('RIS orientation / rad');
    title('g block orientation trajectory');
    legend('\psi_1','\psi_2','\psi_3','Location','best');
    grid on;

    figure(6);
    phi_show_idx = 1:min(8, size(phi_track, 2));
    plot(angle(phi_track(:, phi_show_idx)), 'LineWidth', 1.2);
    xlabel('AO Iteration'); ylabel('Phase angle of diag(\Phi) / rad');
    title('\Phi block variable trajectory (selected RIS elements)');
    legend_str = cell(1, length(phi_show_idx));
    for ii = 1:length(phi_show_idx)
        legend_str{ii} = ['\theta_', num2str(phi_show_idx(ii))];
    end
    legend(legend_str, 'Location', 'best');
    grid on;

    figure(7);
    f_show_idx = 1:min(8, size(f_track, 2));
    plot(abs(f_track(:, f_show_idx)), 'LineWidth', 1.2);
    xlabel('AO Iteration'); ylabel('|f_i|');
    title('Beamforming f variable trajectory (selected coefficients)');
    legend_str = cell(1, length(f_show_idx));
    for ii = 1:length(f_show_idx)
        legend_str{ii} = ['f_', num2str(f_show_idx(ii))];
    end
    legend(legend_str, 'Location', 'best');
    grid on;

    figure(8);
    if ~isempty(loss_track)
        plot(loss_track, '^-', 'LineWidth', 1.2);
        xlabel('AO Iteration'); ylabel('AO objective value');
        title('Overall AO objective after full g/\Phi/f update');
        grid on;
    else
        text(0.1,0.5,'No objective data available.'); axis off;
    end

    figure(9);
    plot(g_delta_history{jj,kk}, 'o-', 'LineWidth', 1.1); hold on;
    plot(phi_delta_history{jj,kk}, 's-', 'LineWidth', 1.1);
    plot(f_delta_history{jj,kk}, 'd-', 'LineWidth', 1.1);
    xlabel('AO Iteration'); ylabel('Block variable update norm');
    title('Update norm of g, \Phi and f blocks');
    legend('||g_t-g_{t-1}||_2','||v_t-v_{t-1}||_2','||f_t-f_{t-1}||_2','Location','best');
    grid on; hold off;
else
    figure(5);
    text(0.1,0.5,'No fixed-N tracking data available.'); axis off;
end

figure(1)
%左轴
yyaxis left;   
plot_j = 3; % 与上面的 for j=3:3 保持一致
if run_single_N
    N_plot = N_total(k_range);
    SNR_plot = SNR_s_dB_final(plot_j, k_range);
    rho_plot = rho_total_final(plot_j, k_range);
else
    valid_idx = find(abs(SNR_s_dB_final(plot_j,:)) > 1e-12 | abs(rho_total_final(plot_j,:)) > 1e-12);
    if isempty(valid_idx), valid_idx = 1:ris_num; end
    N_plot = N_total(valid_idx);
    SNR_plot = SNR_s_dB_final(plot_j, valid_idx);
    rho_plot = rho_total_final(plot_j, valid_idx);
end
plot(N_plot, SNR_plot,"^-");hold on;
% plot(N_total, SNR_s_dB_final(3,:),"^-");hold on;
ylabel('SNR of Sensing');
ylim([0 50]);           % 限制横坐标范围
set(gca, 'YTick', 0:5:60); 
yyaxis right;
plot(rho_plot*0 + N_plot, rho_plot,"^-");hold on;
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