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
Pt_dB = 15; %发射功率为20dBm
Pt = 1e-3*10^(Pt_dB / 10); 
B = 3e6; % 带宽3MHz
sigma_s2_dB = -174 +10*log10(B); %感知场景噪声方差约为-110dBm
sigma_c2_dB = -174 +10*log10(B); %通信场景噪声方差约为-110dBm
sigma_s2 = 1e-3*10^(sigma_s2_dB / 10);
sigma_c2 = 1e-3*10^(sigma_c2_dB / 10);
gamma_0_dB = 10; %通信场景信噪比阈值10dB
gamma_0 = 10^(gamma_0_dB / 10);
%RIS的行反射单元数目
Nx_total = [12]; % RIS每行/列反射单元数，用于多RIS规模速率扫描
N_total = Nx_total.^2;
ris_num = length(Nx_total);
constant_factor = 1e12; % 数值稳定：避免目标/梯度被过度放大导致优化停滞

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
params.phase_objective_scale = constant_factor;
params.outer_objective_scale = 1;
params.secrecy_rate_min = 0; % v2公式10约束阈值（默认无硬门限）
params.secrecy_leakage_weight = 0.5; % 波束对"E端”泄漏抑制权重
params.objective_type = 'v2_formula3_secrecy_rate';
params.mu_v2 = 1.5;
params.grad_eps_phase = 1e-4;
params.manopt_verbosity = 0;
params.manopt_maxiter = 120;
params.manopt_tolgradnorm = 1e-8;
params.manopt_minstepsize = 1e-10;
params.phase_mm_maxiter = 10;        % Phi块内部MM循环次数：theta_k -> beta_k -> theta_{k+1}
params.phase_mm_tol = 1e-4;          % Phi块内部MM循环停止阈值
params.ao_rel_gap_tol = 1e-4;        % AO块间目标差收敛阈值
params.ao_rel_progress_tol = 5e-4;   % AO整轮真实目标相对变化收敛阈值
params.min_ao_iter = 20;             % 至少迭代若干轮，避免过早停在短曲线
params.ao_stall_window = 8;          % 连续若干轮变化都很小才停止
params.max_ao_iter = 500;            % 防止max_iter过大导致未收敛时长时间随机震荡
params.convergence_check_window = 5;    % 额外收敛认证窗口：最后若干轮都足够平稳才标记为已收敛
params.convergence_abs_tol = 1e-4;      % 保密率绝对变化阈值，防止速率接近0时相对误差失真

p_R_begin = [100,100]; %固定点A=(100,100)
location_num = 3; %第一个为A扩展的小区域
%第二个为A扩展的大区域

% 变量个数
num_vars = 5;     
% 定义粒子群大小，变量个数num_vars<= 10时，总数保持在20~50就行;
num_particles = 150; %对不同区域可以适当修改，区域较小的，50即可，较大的200

% 并行设置：use_parallel_sweeps=true时，Alice-Bob距离扫描使用parfor并行。
% 注意：parfor外层并行时，内部PSO不再开UseParallel，避免嵌套并行冲突。
use_parallel_sweeps = true;
parallel_profile = 'threads'; % 优先使用线程池，Processes启动失败时通常可改用threads；若仍失败可设为''
parallel_workers = []; % []表示让MATLAB自动决定worker数量；也可填4、8等
use_pso_parallel = false;
params.pso_use_parallel = use_pso_parallel;

limit = zeros(location_num,ris_num);
SNR_c_dB_final = zeros(location_num,ris_num);
SNR_s_dB_final = zeros(location_num,ris_num);
rho_total_final = zeros(location_num,ris_num);
secrecy_rate_final = nan(location_num,ris_num); % 每个RIS规模最终保密率
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
convergence_status = false(location_num, ris_num); % 每个(j,k)小实验是否通过最终收敛认证
convergence_gap_final = nan(location_num, ris_num); % 每个(j,k)最后窗口内最大相邻保密率变化
convergence_rel_gap_final = nan(location_num, ris_num); % 每个(j,k)最后窗口内最大相对变化
convergence_stop_reason = strings(location_num, ris_num); % stall_window/max_iter/insufficient_history
rho = 0;

% 是否固定一个N做收敛剖析
run_single_N = false;
target_k_track = 1; % run_single_N=true时使用；当前Nx_total(1)=8, N=64

% 是否额外绘制V2论文Fig.2类型曲线：保密速率随Alice-Bob水平距离变化
run_alice_bob_distance_sweep = true;
alice_bob_distance_total = 50:5:80; % Alice-Bob水平距离扫描点，单位m
distance_sweep_k = min(4, ris_num); % 固定一个RIS规模做距离扫描；默认Nx_total(4)=20，即N=400
distance_repeat_num = 1; % 每个Alice-Bob距离重复运行次数，取平均保密率画图
secrecy_rate_vs_distance = nan(size(alice_bob_distance_total));
secrecy_rate_vs_distance_std = nan(size(alice_bob_distance_total));
distance_repeat_rate_samples = nan(length(alice_bob_distance_total), distance_repeat_num);
rho_vs_distance = nan(size(alice_bob_distance_total));
baseline_secrecy_rate_vs_distance = nan(size(alice_bob_distance_total)); % 固定RIS位置/相位/波束的Alice-Bob距离baseline
baseline_rho_vs_distance = nan(size(alice_bob_distance_total));
distance_convergence_status = false(size(alice_bob_distance_total));
distance_convergence_gap_final = nan(size(alice_bob_distance_total));
distance_convergence_rel_gap_final = nan(size(alice_bob_distance_total));
distance_convergence_stop_reason = strings(size(alice_bob_distance_total));
convergence_history_distance = cell(size(alice_bob_distance_total));

% 无人机/6DMA水平移动范围：1=小范围，2=大范围，3=固定位置。
% 若要画完整“移动范围-保密性能”对比图，可改成 j_range = 1:3。
j_range = 1:3;

if use_parallel_sweeps
    use_parallel_sweeps = ensure_parallel_pool(parallel_workers);
end

for j=j_range %一个一个跑效果更好
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
        max_iter = params.max_ao_iter; % 最大迭代次数上限，实际停止由AO收敛条件决定
        % 初始化v
        theta = linspace(0,pi,N);
        v_init = exp(1j * theta);
        v = v_init';
        vetor_total=[];

        init_p_R = x_opt(1:2).';
        init_psi = x_opt(3:5).';
        [~,~, hc_init, ht_init, ~, ~, ~, ~, ~, ~, ~] = bulid_H(init_p_R, init_psi, v, params);
        f_current = solve_beamforming_vector(hc_init, ht_init, params);
        params.fixed_beamformer = f_current;
        [~, initial_secrecy_rate, ~, ~, ~] = secrecy_rate_objective(hc_init, ht_init, params);
        params = rmfield(params, 'fixed_beamformer');

        f_val_2_pre = initial_secrecy_rate;
        stall_count = 0;
        stop_reason = "max_iter";
        f_curve = [];
        secrecy_rate_curve = initial_secrecy_rate; % 第0次AO迭代：初始g/Phi/f对应的真实保密率
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
            'UseParallel', params.pso_use_parallel, ...
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
        
            % 固定 p_R, psi，使用V2式(10)的相位MM子问题优化v
            params.phase_surrogate_mode = true;
            [v, f_val_2, info] = rieman_grad_fun(p_R, psi, x_opt, v, params);
            params = rmfield(params, 'phase_surrogate_mode');
            g_curve = [g_curve; f_val_1];
            phi_curve = [phi_curve; f_val_2];
            phi_value_curve = [phi_value_curve; v(:).'];
            phi_delta_curve = [phi_delta_curve; norm(v - v_prev_for_delta, 2)];

            % f块：固定当前g和Phi后，重新构造/更新波束，并记录更新后的真实目标
            [rho_iter,~, hc_tmp, ht_tmp, ~, Uc_tmp, a_BU_tmp, Ut_tmp, a_BT_tmp, ~, ~] = bulid_H(p_R, psi, v, params);
            f_tmp = solve_beamforming_vector(hc_tmp, ht_tmp, params);
            params.fixed_beamformer = f_tmp;
            [f_obj_raw,secrecy_rate_iter,~,~,~] = secrecy_rate_objective(hc_tmp, ht_tmp, params);
            f_obj = params.constant_factor * f_obj_raw;

            f_curve_beam = [f_curve_beam; f_obj];
            loss_curve = [loss_curve; f_obj];
            f_curve = [f_curve; secrecy_rate_iter];
            secrecy_rate_curve = [secrecy_rate_curve; secrecy_rate_iter];
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
                num2str(f_val_1), ', Phi更新后函数值: ', num2str(f_val_2), ', f更新后函数值: ', num2str(f_obj), ...
                ', 保密率: ', num2str(secrecy_rate_iter), ',信道相关性为：',num2str(rho_iter)]);
            rho = rho_iter;
        
             % 保存当前值
            p_R_pre = p_R;
            psi_final = psi;
            v_final = v; 
        
            % 可选：检查是否收敛或达到停止条件
            %rel_gap = abs(f_val_1 - f_obj) / max(abs(f_val_1), 1e-12);
            if isfinite(f_val_2_pre)
                rel_progress = abs(secrecy_rate_iter - f_val_2_pre) / max(abs(f_val_2_pre), 1e-12);
            else
                rel_progress = inf;
            end
            if iter >= params.min_ao_iter && rel_progress < params.ao_rel_progress_tol
                stall_count = stall_count + 1;
            else
                stall_count = 0;
            end
            if stall_count >= params.ao_stall_window
                stop_reason = "stall_window";
                break;
            end
            f_val_2_pre = secrecy_rate_iter;

        end
        if isempty(f_curve)
            f_curve = f_val_2_pre;
        end
        convergence_history{j, k} = secrecy_rate_curve;
        g_obj_history{j, k} = g_curve;
        phi_obj_history{j, k} = phi_curve;
        f_obj_history{j, k} = f_curve_beam;
        loss_history{j, k} = loss_curve;
        total_obj_history{j, k} = loss_curve;
        if ~isempty(secrecy_rate_curve)
            secrecy_rate_final(j,k) = secrecy_rate_curve(end);
        end
        g_value_history{j, k} = g_value_curve;
        phi_value_history{j, k} = phi_value_curve;
        f_value_history{j, k} = f_value_curve;
        g_delta_history{j, k} = g_delta_curve;
        phi_delta_history{j, k} = phi_delta_curve;
        f_delta_history{j, k} = f_delta_curve;
        [is_converged_tmp, abs_gap_tmp, rel_gap_tmp] = certify_convergence(secrecy_rate_curve, params);
        convergence_status(j,k) = is_converged_tmp;
        convergence_gap_final(j,k) = abs_gap_tmp;
        convergence_rel_gap_final(j,k) = rel_gap_tmp;
        if ~is_converged_tmp && stop_reason == "stall_window"
            stop_reason = "failed_final_certificate";
        end
        convergence_stop_reason(j,k) = stop_reason;
        disp(['收敛认证 j=', num2str(j), ', k=', num2str(k), ', N=', num2str(N), ...
            ': converged=', mat2str(is_converged_tmp), ', final_abs_gap=', num2str(abs_gap_tmp), ...
            ', final_rel_gap=', num2str(rel_gap_tmp), ', stop_reason=', char(stop_reason)]);

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

%% V2论文Fig.2类型图：保密速率随Alice-Bob水平距离变化
% 说明：原始RIS规模扫描只改变N；这里固定一个RIS规模，逐个改变Bob/UE的x坐标，
% 从而得到 R_s versus horizontal distance between Alice and Bob。
if run_alice_bob_distance_sweep
    distance_plot_j = j;
    distance_Nx = Nx_total(distance_sweep_k);
    rng(45);
    if distance_plot_j == 1
        distance_x_lb = [50, 50, 0, 0, 0];
        distance_x_ub = [100, 100, 2*pi, 2*pi, 2*pi];
        distance_x_init = distance_x_lb + (distance_x_ub - distance_x_lb) .* rand(1, num_vars);
        distance_x_init(1:2) = p_R_begin;
        distance_x_init(3:end) = [3.74, 2, 5.98];
    elseif distance_plot_j == 2
        distance_x_lb = [0, 0, 0, 0, 0];
        distance_x_ub = [100, 100, 2*pi, 2*pi, 2*pi];
        distance_x_init = distance_x_lb + (distance_x_ub - distance_x_lb) .* rand(1, num_vars);
        distance_x_init(1:2) = p_R_begin;
        distance_x_init(3:end) = [3.74, 2, 5.98];
    else
        distance_x_lb = [100, 100, 0, 0, 0];
        distance_x_ub = [100, 100, 2*pi, 2*pi, 2*pi];
        distance_x_init = distance_x_lb + (distance_x_ub - distance_x_lb) .* rand(1, num_vars);
        distance_x_init(1:2) = p_R_begin;
        distance_x_init(3:end) = [3.74, 2, 5.98];
    end

    % 固定RIS位置、固定相位、固定波束的baseline：
    % 只改变Bob/UE的水平距离，不重新优化g/Phi/f，用于观察纯距离变化趋势。
    params_baseline_distance = params;
    params_baseline_distance.Nx = distance_Nx;
    params_baseline_distance.Ny = distance_Nx;
    params_baseline_distance.N = distance_Nx^2;
    theta_baseline_distance = linspace(0, pi, params_baseline_distance.N);
    v_baseline_distance = exp(1j * theta_baseline_distance).';
    p_R_baseline_distance = distance_x_init(1:2).';
    psi_baseline_distance = distance_x_init(3:5).';
    params_baseline_distance.p_U = [alice_bob_distance_total(1), 0, 0].';
    [~,~, hc_baseline_ref, ht_baseline_ref, ~, ~, ~, ~, ~, ~, ~] = ...
        bulid_H(p_R_baseline_distance, psi_baseline_distance, v_baseline_distance, params_baseline_distance);
    f_baseline_distance = solve_beamforming_vector(hc_baseline_ref, ht_baseline_ref, params_baseline_distance);

    if use_parallel_sweeps
        parfor d_idx = 1:length(alice_bob_distance_total)
            [secrecy_rate_vs_distance(d_idx), secrecy_rate_vs_distance_std(d_idx), rho_vs_distance(d_idx), ...
                convergence_history_distance{d_idx}, distance_convergence_status(d_idx), ...
                distance_convergence_gap_final(d_idx), distance_convergence_rel_gap_final(d_idx), ...
                distance_convergence_stop_reason(d_idx), baseline_secrecy_rate_vs_distance(d_idx), ...
                baseline_rho_vs_distance(d_idx), distance_repeat_rate_samples(d_idx, :)] = ...
                run_distance_point_repeated(alice_bob_distance_total(d_idx), params, params_baseline_distance, ...
                p_R_baseline_distance, psi_baseline_distance, v_baseline_distance, f_baseline_distance, ...
                distance_x_lb, distance_x_ub, distance_x_init, distance_Nx, num_particles, num_vars, distance_repeat_num);
        end
    else
        for d_idx = 1:length(alice_bob_distance_total)
            [secrecy_rate_vs_distance(d_idx), secrecy_rate_vs_distance_std(d_idx), rho_vs_distance(d_idx), ...
                convergence_history_distance{d_idx}, distance_convergence_status(d_idx), ...
                distance_convergence_gap_final(d_idx), distance_convergence_rel_gap_final(d_idx), ...
                distance_convergence_stop_reason(d_idx), baseline_secrecy_rate_vs_distance(d_idx), ...
                baseline_rho_vs_distance(d_idx), distance_repeat_rate_samples(d_idx, :)] = ...
                run_distance_point_repeated(alice_bob_distance_total(d_idx), params, params_baseline_distance, ...
                p_R_baseline_distance, psi_baseline_distance, v_baseline_distance, f_baseline_distance, ...
                distance_x_lb, distance_x_ub, distance_x_init, distance_Nx, num_particles, num_vars, distance_repeat_num);
            disp(['Alice-Bob水平距离 ', num2str(alice_bob_distance_total(d_idx)), ' m, RIS规模N=', ...
                num2str(distance_Nx^2), ', 最终保密率: ', num2str(secrecy_rate_vs_distance(d_idx)), ...
                ', 信道相关性: ', num2str(rho_vs_distance(d_idx)), ...
                ', baseline保密率: ', num2str(baseline_secrecy_rate_vs_distance(d_idx)), ...
                ', 重复均值标准差: ', num2str(secrecy_rate_vs_distance_std(d_idx)), ...
                ', converged=', mat2str(distance_convergence_status(d_idx)), ...
                ', final_rel_gap=', num2str(distance_convergence_rel_gap_final(d_idx))]);
        end
    end
    params.p_U = p_U;
end

% 收敛折线图（示例：当前运行的 j=3，k=4）
figure(2);
nonempty_mask = ~cellfun(@isempty, convergence_history);
[j_idx, k_idx] = find(nonempty_mask);
if ~isempty(j_idx)
    % 自动选取最后一个有数据的(j,k)，避免硬编码索引没有数据
    target_j = j_idx(end);
    target_k = k_idx(end);
    iter_axis = 0:(length(convergence_history{target_j, target_k}) - 1);
    plot(iter_axis, convergence_history{target_j, target_k}, 'o-', 'LineWidth', 1.2);
    xlabel('AO Iteration');
    ylabel('Secrecy Rate (bit/s/Hz)');
    title(['Convergence Curve (j=', num2str(target_j), ', k=', num2str(target_k), ')']);
    grid on;
    % 图2-2：相邻迭代目标变化（观察是否进入平台期）
    if length(convergence_history{target_j, target_k}) >= 2
        figure(3);
        delta_obj = abs(diff(convergence_history{target_j, target_k}));
        gap_axis = 1:(length(convergence_history{target_j, target_k}) - 1);
        semilogy(gap_axis, delta_obj, 's-', 'LineWidth', 1.1);
        xlabel('AO Iteration');
        ylabel('|R_s(t) - R_s(t-1)| (log scale)');
        title(['Convergence gap (j=', num2str(target_j), ', k=', num2str(target_k), ')']);
        grid on;
    end

%     % 图2-3：同一位置下，不同RIS规模的收敛曲线对比（最多4条）。
%     % run_single_N=true时只会有一个k有数据，因此只画实际存在的数据，避免空白图。
%     figure(4);
%     hold on;
%     available_k = find(~cellfun(@isempty, convergence_history(target_j, :)));
%     k_show = available_k;
%     legend_str = {};
%     for kk = k_show
%         plot(convergence_history{target_j, kk}, 'LineWidth', 1.1);
%         legend_str{end+1} = ['N=', num2str((Nx_total(kk))^2)]; %#ok<SAGROW>
%     end
%     xlabel('AO Iteration');
%     ylabel('Secrecy Rate (bit/s/Hz)');
%     title(['Convergence comparison across available RIS sizes (j=', num2str(target_j), ')']);
%     grid on;
%     if ~isempty(legend_str)
%         legend(legend_str, 'Location', 'best');
%     else
%         text(0.1,0.5,'No RIS-size comparison data available.'); axis off;
%     end
%     hold off;
% else
%     text(0.1, 0.5, 'No convergence data generated in this run.');
%     axis off;
% end
% 
% % 画一个代表性RIS规模下三个块变量自身的变化；多RIS扫描时默认取最后一个有数据的规模
% if ~isempty(j_idx)
%     jj = j_idx(end);
%     kk = k_idx(end);
% elseif run_single_N
%     jj = j;
%     kk = k_range(1);
% else
%     jj = [];
%     kk = [];
% end
% if ~isempty(jj) && ~isempty(g_value_history{jj,kk})
%     g_track = g_value_history{jj,kk};
%     phi_track = phi_value_history{jj,kk};
%     f_track = f_value_history{jj,kk};
%     loss_track = total_obj_history{jj,kk};
% 
%     figure(5);
%     subplot(2,1,1);
%     plot(g_track(:,1:2), 'LineWidth', 1.2);
%     xlabel('AO Iteration'); ylabel('RIS position / m');
%     title(['g block position trajectory (N=', num2str(Nx_total(kk)^2), ')']);
%     legend('p_R_x','p_R_y','Location','best');
%     grid on;
% 
%     subplot(2,1,2);
%     plot(g_track(:,3:5), 'LineWidth', 1.2);
%     xlabel('AO Iteration'); ylabel('RIS orientation / rad');
%     title('g block orientation trajectory');
%     legend('\psi_1','\psi_2','\psi_3','Location','best');
%     grid on;
% 
%     figure(6);
%     phi_show_idx = 1:min(8, size(phi_track, 2));
%     plot(angle(phi_track(:, phi_show_idx)), 'LineWidth', 1.2);
%     xlabel('AO Iteration'); ylabel('Phase angle of diag(\Phi) / rad');
%     title('\Phi block variable trajectory (selected RIS elements)');
%     legend_str = cell(1, length(phi_show_idx));
%     for ii = 1:length(phi_show_idx)
%         legend_str{ii} = ['\theta_', num2str(phi_show_idx(ii))];
%     end
%     legend(legend_str, 'Location', 'best');
%     grid on;
% 
%     figure(7);
%     f_show_idx = 1:min(8, size(f_track, 2));
%     plot(abs(f_track(:, f_show_idx)), 'LineWidth', 1.2);
%     xlabel('AO Iteration'); ylabel('|f_i|');
%     title('Beamforming f variable trajectory (selected coefficients)');
%     legend_str = cell(1, length(f_show_idx));
%     for ii = 1:length(f_show_idx)
%         legend_str{ii} = ['f_', num2str(f_show_idx(ii))];
%     end
%     legend(legend_str, 'Location', 'best');
%     grid on;
% 
%     figure(8);
%     if ~isempty(loss_track)
%         plot(loss_track, '^-', 'LineWidth', 1.2);
%         xlabel('AO Iteration'); ylabel('AO objective value');
%         title('Overall AO objective after full g/\Phi/f update');
%         grid on;
%     else
%         text(0.1,0.5,'No objective data available.'); axis off;
%     end
% 
%     figure(9);
%     phi_delta_per_element = phi_delta_history{jj,kk} / sqrt(max(size(phi_track, 2), 1));
%     plot(g_delta_history{jj,kk}, 'o-', 'LineWidth', 1.1); hold on;
%     plot(phi_delta_per_element, 's-', 'LineWidth', 1.1);
%     plot(f_delta_history{jj,kk}, 'd-', 'LineWidth', 1.1);
%     xlabel('AO Iteration'); ylabel('Block variable update norm');
%     title('Update norm of g, average \Phi element and f blocks');
%     legend('||g_t-g_{t-1}||_2','||v_t-v_{t-1}||_2 / sqrt(N)','||f_t-f_{t-1}||_2','Location','best');
%     grid on; hold off;
% else
%     figure(5);
%     text(0.1,0.5,'No block-tracking data available.'); axis off;
% end
% 
% % 图10：RIS规模-最终保密率曲线，用来直观看RIS数量变化后的模型效果
% figure(10);
% rate_plot_j = j;
% valid_rate_idx = find(~isnan(secrecy_rate_final(rate_plot_j,:)));
% if ~isempty(valid_rate_idx)
%     plot(Nx_total(valid_rate_idx), secrecy_rate_final(rate_plot_j, valid_rate_idx), 'o-', 'LineWidth', 1.4);
%     xlabel('RIS elements per side N_x=N_y');
%     ylabel('Secrecy Rate (bit/s/Hz)');
%     title(['Final secrecy rate versus RIS size (N=N_x^2, j=', num2str(rate_plot_j), ')']);
%     grid on;
% else
%     text(0.1,0.5,'No secrecy-rate data available.'); axis off;
% end
% 
% 
% 
% % 图11：总RIS反射单元数N=Nx*Ny与最终保密率的关系（论文主图之一）
% figure(11);
% valid_rate_idx = find(~isnan(secrecy_rate_final(rate_plot_j,:)));
% if ~isempty(valid_rate_idx)
%     plot(N_total(valid_rate_idx), secrecy_rate_final(rate_plot_j, valid_rate_idx), 's-', 'LineWidth', 1.4);
%     xlabel('Total number of RIS reflecting elements N=N_xN_y');
%     ylabel('Secrecy Rate (bit/s/Hz)');
%     title(['Final secrecy rate versus total RIS elements (j=', num2str(rate_plot_j), ')']);
%     grid on;
% else
%     text(0.1,0.5,'No secrecy-rate data available.'); axis off;
% end
% 
% % 图12：不同RIS规模下的归一化收敛增益，突出每条曲线自身的提升过程
% figure(12);
% hold on;
% legend_str = {};
% if exist('target_j', 'var') && ~isempty(target_j)
%     available_k_norm = find(~cellfun(@isempty, convergence_history(target_j, :)));
% else
%     available_k_norm = valid_rate_idx;
%     target_j = rate_plot_j;
% end
% for kk = available_k_norm
%     curve_tmp = convergence_history{target_j, kk};
%     if ~isempty(curve_tmp)
%         plot(curve_tmp - curve_tmp(1), 'LineWidth', 1.1);
%         legend_str{end+1} = ['N=', num2str(N_total(kk))]; %#ok<SAGROW>
%     end
% end
% xlabel('AO Iteration');
% ylabel('Secrecy-rate gain R_s(t)-R_s(1) (bit/s/Hz)');
% title(['Normalized convergence gain across RIS sizes (j=', num2str(target_j), ')']);
% % 注意：图12展示的是累计增益，不是收敛判据；真正的收敛认证见图16。
% grid on;
% if ~isempty(legend_str)
%     legend(legend_str, 'Location', 'best');
% else
%     text(0.1,0.5,'No convergence data available.'); axis off;
% end
% hold off;
% 
% % 图13：最终保密率柱状图，适合论文中直观比较不同RIS规模
% figure(13);
% if ~isempty(valid_rate_idx)
%     bar(Nx_total(valid_rate_idx), secrecy_rate_final(rate_plot_j, valid_rate_idx));
%     xlabel('RIS elements per side N_x=N_y');
%     ylabel('Final Secrecy Rate (bit/s/Hz)');
%     title(['Bar chart of final secrecy rate across RIS sizes (j=', num2str(rate_plot_j), ')']);
%     grid on;
% else
%     text(0.1,0.5,'No secrecy-rate data available.'); axis off;
% end
% 
% % 图14：保密率、感知SNR和相关系数rho的联合趋势，展示隐私目标与ISAC诊断量
% figure(14);
% if ~isempty(valid_rate_idx)
%     yyaxis left;
%     plot(N_total(valid_rate_idx), secrecy_rate_final(rate_plot_j, valid_rate_idx), 'o-', 'LineWidth', 1.3); hold on;
%     plot(N_total(valid_rate_idx), SNR_s_dB_final(rate_plot_j, valid_rate_idx), '^-', 'LineWidth', 1.3);
%     ylabel('Rate / Sensing SNR(dB)');
%     yyaxis right;
%     plot(N_total(valid_rate_idx), rho_total_final(rate_plot_j, valid_rate_idx), 's--', 'LineWidth', 1.3);
%     ylabel('\rho');
%     xlabel('Total number of RIS reflecting elements N=N_xN_y');
%     title(['Final metrics versus RIS size (j=', num2str(rate_plot_j), ')']);
%     legend('Secrecy Rate','Sensing SNR','\rho','Location','best');
%     grid on;
% else
%     text(0.1,0.5,'No final metric data available.'); axis off;
% end
% 
% % 图15：对标V2 Fig.2，Alice-Bob水平距离与最终保密速率关系
% figure(15);
% valid_distance_idx = find(~isnan(secrecy_rate_vs_distance));
% if ~isempty(valid_distance_idx)
%     plot(alice_bob_distance_total(valid_distance_idx), secrecy_rate_vs_distance(valid_distance_idx), ...
%         'o-', 'LineWidth', 1.4); hold on;
%     legend_items = {'AO optimized secrecy rate (best of repeats)'};
%     valid_baseline_distance_idx = find(~isnan(baseline_secrecy_rate_vs_distance));
%     if ~isempty(valid_baseline_distance_idx)
%         plot(alice_bob_distance_total(valid_baseline_distance_idx), ...
%             baseline_secrecy_rate_vs_distance(valid_baseline_distance_idx), 's--', 'LineWidth', 1.3);
%         legend_items{end+1} = 'Fixed g/Phi/f baseline'; %#ok<SAGROW>
%     end
%     unconverged_distance_idx = valid_distance_idx(~distance_convergence_status(valid_distance_idx));
%     if ~isempty(unconverged_distance_idx)
%         plot(alice_bob_distance_total(unconverged_distance_idx), secrecy_rate_vs_distance(unconverged_distance_idx), ...
%             'rx', 'MarkerSize', 9, 'LineWidth', 1.4);
%         legend_items{end+1} = 'Not certified converged'; %#ok<SAGROW>
%     end
%     legend(legend_items, 'Location', 'best');
%     xlabel('Horizontal distance between Alice and Bob (m)');
%     ylabel('Secrecy Rate (bit/s/Hz)');
%     title(['Secrecy rate versus Alice-Bob horizontal distance (N=', num2str(Nx_total(distance_sweep_k)^2), ')']);
%     grid on; hold off;
% else
%     text(0.1,0.5,'No Alice-Bob distance sweep data available.'); axis off;
% end
% 
% % 图16：每个RIS规模小实验的最终收敛认证，低于阈值且状态为true才建议用于论文结果
% figure(16);
% valid_conv_idx = find(~isnan(convergence_rel_gap_final(rate_plot_j,:)));
% if ~isempty(valid_conv_idx)
%     semilogy(N_total(valid_conv_idx), convergence_rel_gap_final(rate_plot_j, valid_conv_idx), 'o-', 'LineWidth', 1.3); hold on;
%     yline(params.ao_rel_progress_tol, 'r--', 'Convergence threshold');
%     not_certified_idx = valid_conv_idx(~convergence_status(rate_plot_j, valid_conv_idx));
%     if ~isempty(not_certified_idx)
%         semilogy(N_total(not_certified_idx), convergence_rel_gap_final(rate_plot_j, not_certified_idx), ...
%             'rx', 'MarkerSize', 9, 'LineWidth', 1.4);
%         legend('Final relative gap', 'Threshold', 'Not certified converged', 'Location', 'best');
%     else
%         legend('Final relative gap', 'Threshold', 'Location', 'best');
%     end
%     xlabel('Total number of RIS reflecting elements N=N_xN_y');
%     ylabel('Certified final relative gap');
%     title(['Convergence certificate across RIS sizes (j=', num2str(rate_plot_j), ')']);
%     grid on; hold off;
% else
%     text(0.1,0.5,'No convergence certificate data available.'); axis off;
% end
% 
% % 图17：Alice-Bob距离扫描的最终收敛认证
% figure(17);
% valid_distance_conv_idx = find(~isnan(distance_convergence_rel_gap_final));
% if ~isempty(valid_distance_conv_idx)
%     semilogy(alice_bob_distance_total(valid_distance_conv_idx), distance_convergence_rel_gap_final(valid_distance_conv_idx), ...
%         's-', 'LineWidth', 1.3); hold on;
%     yline(params.ao_rel_progress_tol, 'r--', 'Convergence threshold');
%     not_certified_distance_idx = valid_distance_conv_idx(~distance_convergence_status(valid_distance_conv_idx));
%     if ~isempty(not_certified_distance_idx)
%         semilogy(alice_bob_distance_total(not_certified_distance_idx), ...
%             distance_convergence_rel_gap_final(not_certified_distance_idx), 'rx', 'MarkerSize', 9, 'LineWidth', 1.4);
%         legend('Final relative gap', 'Threshold', 'Not certified converged', 'Location', 'best');
%     else
%         legend('Final relative gap', 'Threshold', 'Location', 'best');
%     end
%     xlabel('Horizontal distance between Alice and Bob (m)');
%     ylabel('Certified final relative gap');
%     title('Convergence certificate for Alice-Bob distance sweep');
%     grid on; hold off;
% else
%     text(0.1,0.5,'No distance convergence certificate data available.'); axis off;
% end
% 
% % 图18：无人机/6DMA移动范围与保密性能对比
% figure(18);
% mobility_plot_order = [3, 1, 2]; % 按“固定位置 -> 小范围 -> 大范围”的活动范围递增顺序画图
% mobility_labels_ordered = {'固定位置', '小范围移动', '大范围移动'};
% mobility_rate_raw = nan(1, length(mobility_plot_order));
% for idx_tmp = 1:length(mobility_plot_order)
%     jj_tmp = mobility_plot_order(idx_tmp);
%     rate_candidates = secrecy_rate_final(jj_tmp, :);
%     valid_candidates = find(~isnan(rate_candidates));
%     if ~isempty(valid_candidates)
%         mobility_rate_raw(idx_tmp) = rate_candidates(valid_candidates(end));
%     end
% end
% if any(~isnan(mobility_rate_raw))
%     mobility_rate_envelope = mobility_rate_raw;
%     for idx_tmp = 2:length(mobility_rate_envelope)
%         mobility_rate_envelope(idx_tmp) = max(mobility_rate_envelope(1:idx_tmp), [], 'omitnan');
%     end
%     x_mobility = categorical(mobility_labels_ordered, mobility_labels_ordered, 'Ordinal', true);
%     bar(x_mobility, [mobility_rate_raw(:), mobility_rate_envelope(:)]);
%     ylabel('Final Secrecy Rate (bit/s/Hz)');
%     title('Secrecy performance versus UAV/6DMA movement range');
%     legend('AO result', 'Best feasible reuse/envelope', 'Location', 'best');
%     grid on;
% else
%     text(0.1,0.5,'No movement-range comparison data available. Set j_range = 1:3 to compare ranges.'); axis off;
% end
% 
% % 图18：无人机/6DMA移动范围与保密性能对比
% figure(18);
% mobility_labels = {'小范围移动', '大范围移动', '固定位置'};
% valid_mobility_idx = find(any(~isnan(secrecy_rate_final), 2)).';
% if ~isempty(valid_mobility_idx)
%     mobility_rate = nan(size(valid_mobility_idx));
%     for idx_tmp = 1:length(valid_mobility_idx)
%         jj_tmp = valid_mobility_idx(idx_tmp);
%         rate_candidates = secrecy_rate_final(jj_tmp, :);
%         valid_candidates = find(~isnan(rate_candidates));
%         if ~isempty(valid_candidates)
%             mobility_rate(idx_tmp) = rate_candidates(valid_candidates(end));
%         end
%     end
%     bar(categorical(mobility_labels(valid_mobility_idx)), mobility_rate);
%     ylabel('Final Secrecy Rate (bit/s/Hz)');
%     title('Secrecy performance versus UAV/6DMA movement range');
%     grid on;
% else
%     text(0.1,0.5,'No movement-range comparison data available. Set j_range = 1:3 to compare ranges.'); axis off;
end

figure(1)
%左轴
yyaxis left;   
plot_j = j; % 与上面的 for j=3:3 保持一致
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
legend(['j=', num2str(plot_j), ', sensing SNR'], ['j=', num2str(plot_j), ', \rho'], 'Location', 'best');
xlabel('Number of RIS reflection units');
title('Sensing SNR and channel correlation versus RIS size');
  
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
%% 自动保存所有输出图像：每次运行按时间戳新建文件夹，文件名使用中文说明
save_fig_timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
save_fig_dir = fullfile(root_dir, '输出图像', save_fig_timestamp);
if ~exist(save_fig_dir, 'dir')
    mkdir(save_fig_dir);
end
figure_name_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
figure_name_map(1) = '感知SNR与信道相关性随RIS规模变化';
figure_name_map(2) = '单个RIS规模下保密速率收敛曲线';
figure_name_map(3) = '单个RIS规模下保密速率收敛间隔';
% figure_name_map(4) = '不同RIS规模下保密速率收敛对比';
% figure_name_map(5) = 'g块位置和姿态变量轨迹';
% figure_name_map(6) = 'Phi相位变量轨迹';
% figure_name_map(7) = '波束向量f变量轨迹';
% figure_name_map(8) = 'AO整体目标函数变化';
% figure_name_map(9) = 'g_Phi_f三块变量更新范数';
% figure_name_map(10) = '每边RIS单元数与最终保密速率';
% figure_name_map(11) = '总RIS反射单元数与最终保密速率';
% figure_name_map(12) = '不同RIS规模下归一化收敛增益';
% figure_name_map(13) = '不同RIS规模最终保密速率柱状图';
% figure_name_map(14) = '保密速率感知SNR和rho联合趋势';
% figure_name_map(15) = 'Alice_Bob水平距离与保密速率关系';
% figure_name_map(16) = '不同RIS规模收敛认证';
% figure_name_map(17) = 'Alice_Bob距离扫描收敛认证';
% figure_name_map(18) = '无人机移动范围与保密性能对比';
fig_handles = findall(0, 'Type', 'figure');
[~, fig_order] = sort(arrayfun(@(h) h.Number, fig_handles));
fig_handles = fig_handles(fig_order);
for fig_idx = 1:length(fig_handles)
    fig_handle = fig_handles(fig_idx);
    fig_number = fig_handle.Number;
    if isKey(figure_name_map, fig_number)
        fig_name = figure_name_map(fig_number);
    else
        fig_name = ['图', num2str(fig_number)];
    end
    file_prefix = sprintf('图%02d_%s', fig_number, fig_name);
    savefig(fig_handle, fullfile(save_fig_dir, [file_prefix, '.fig']));
    saveas(fig_handle, fullfile(save_fig_dir, [file_prefix, '.png']));
end
% save(fullfile(save_fig_dir, '本次运行数据.mat'), 'Nx_total', 'N_total', 'secrecy_rate_final', ...
%     'SNR_s_dB_final', 'SNR_c_dB_final', 'rho_total_final', 'convergence_history', ...
%     'g_value_history', 'phi_value_history', 'f_value_history', 'alice_bob_distance_total', ...
%     'secrecy_rate_vs_distance', 'secrecy_rate_vs_distance_std', 'distance_repeat_rate_samples', ...
%     'rho_vs_distance', 'baseline_secrecy_rate_vs_distance', ...
%     'baseline_rho_vs_distance', 'convergence_history_distance', ...
%     'convergence_status', 'convergence_gap_final', 'convergence_rel_gap_final', ...
%     'convergence_stop_reason', 'distance_convergence_status', 'distance_convergence_gap_final', ...
%     'distance_convergence_rel_gap_final', 'distance_convergence_stop_reason', ...
%     'j_range', 'use_parallel_sweeps', 'use_pso_parallel', 'distance_repeat_num');
% disp(['所有图像和运行数据已保存到: ', save_fig_dir]);


function use_parallel = ensure_parallel_pool(parallel_workers)
%ENSURE_PARALLEL_POOL 若可用则启动并行池，供parfor距离扫描使用。
use_parallel = false;
if isempty(ver('parallel'))
    warning('未检测到 Parallel Computing Toolbox，自动退回普通for循环。');
    return;
end
try
    pool = gcp('nocreate');
    if isempty(pool)
        if isempty(parallel_workers)
            parpool;
        else
            parpool(parallel_workers);
        end
    end
    use_parallel = true;
catch ME
    warning('并行池启动失败，自动退回普通for循环：%s', ME.message);
end
end

function [mean_rate, std_rate, mean_rho, convergence_curves, all_converged, max_abs_gap, max_rel_gap, stop_reason, baseline_rate, baseline_rho, rate_samples] = ...
    run_distance_point_repeated(distance_value, params, params_baseline_distance, p_R_baseline_distance, ...
    psi_baseline_distance, v_baseline_distance, f_baseline_distance, distance_x_lb, distance_x_ub, ...
    distance_x_init, distance_Nx, num_particles, num_vars, repeat_num)
%RUN_DISTANCE_POINT_REPEATED 同一Alice-Bob距离重复运行多次，返回平均保密率。
rate_samples = nan(1, repeat_num);
rho_samples = nan(1, repeat_num);
abs_gap_samples = nan(1, repeat_num);
rel_gap_samples = nan(1, repeat_num);
converged_samples = false(1, repeat_num);
stop_reason_samples = strings(1, repeat_num);
convergence_curves = cell(1, repeat_num);
baseline_rate = nan;
baseline_rho = nan;

for repeat_idx = 1:repeat_num
    rng(1000 + round(distance_value) * 10 + repeat_idx);
    [rate_samples(repeat_idx), rho_samples(repeat_idx), convergence_curves{repeat_idx}, ...
        converged_samples(repeat_idx), abs_gap_samples(repeat_idx), rel_gap_samples(repeat_idx), ...
        stop_reason_samples(repeat_idx), baseline_rate_tmp, baseline_rho_tmp] = ...
        run_distance_point(distance_value, params, params_baseline_distance, p_R_baseline_distance, ...
        psi_baseline_distance, v_baseline_distance, f_baseline_distance, distance_x_lb, distance_x_ub, ...
        distance_x_init, distance_Nx, num_particles, num_vars);
    if repeat_idx == 1
        baseline_rate = baseline_rate_tmp;
        baseline_rho = baseline_rho_tmp;
    end
end

mean_rate = mean(rate_samples, 'omitnan');
std_rate = std(rate_samples, 0, 'omitnan');
mean_rho = mean(rho_samples, 'omitnan');
all_converged = all(converged_samples);
max_abs_gap = max(abs_gap_samples);
max_rel_gap = max(rel_gap_samples);
if all_converged
    stop_reason = "all_repeats_converged";
else
    stop_reason = strjoin(unique(stop_reason_samples), ',');
end
end

function [secrecy_rate_opt, rho_opt, convergence_curve, is_converged, abs_gap, rel_gap, stop_reason, baseline_rate, baseline_rho] = ...
    run_distance_point(distance_value, params, params_baseline_distance, p_R_baseline_distance, ...
    psi_baseline_distance, v_baseline_distance, f_baseline_distance, distance_x_lb, distance_x_ub, ...
    distance_x_init, distance_Nx, num_particles, num_vars)
%DISTANCE_POINT 固定一个Alice-Bob距离，计算baseline并运行一次AO距离扫描。
params_baseline_local = params_baseline_distance;
params_baseline_local.p_U = [distance_value, 0, 0].';
[baseline_rho,~, hc_baseline_tmp, ht_baseline_tmp, ~, ~, ~, ~, ~, ~, ~] = ...
    bulid_H(p_R_baseline_distance, psi_baseline_distance, v_baseline_distance, params_baseline_local);
params_baseline_local.fixed_beamformer = f_baseline_distance;
[~, baseline_rate, ~, ~, ~] = secrecy_rate_objective(hc_baseline_tmp, ht_baseline_tmp, params_baseline_local);

params_distance = params;
params_distance.p_U = [distance_value, 0, 0].';
params_distance.pso_use_parallel = false; % parfor外层并行时避免PSO嵌套并行
[secrecy_rate_opt, rho_opt, convergence_curve, is_converged, abs_gap, rel_gap, stop_reason] = ...
    run_distance_sweep_ao(params_distance, distance_x_lb, distance_x_ub, distance_x_init, ...
    distance_Nx, num_particles, num_vars);
end

function [final_secrecy_rate, final_rho, secrecy_rate_curve, is_converged, final_abs_gap, final_rel_gap, stop_reason] = run_distance_sweep_ao(params_distance, x_lb, x_ub, x_init, Nx, num_particles, num_vars)
%RUN_DISTANCE_SWEEP_AO 固定RIS规模和Alice-Bob距离时，复用主AO流程得到最终保密率。
% 这个局部函数只记录V2 Fig.2所需的最终保密率，不改变主RIS规模扫描结果。
Ny = Nx;
N = Nx * Ny;
params_distance.Nx = Nx;
params_distance.Ny = Ny;
params_distance.N = N;

max_iter = params_distance.max_ao_iter;
theta = linspace(0, pi, N);
v = exp(1j * theta).';
x_opt_distance = x_init;

init_p_R = x_opt_distance(1:2).';
init_psi = x_opt_distance(3:5).';
[~,~, hc_init, ht_init, ~, ~, ~, ~, ~, ~, ~] = bulid_H(init_p_R, init_psi, v, params_distance);
f_current = solve_beamforming_vector(hc_init, ht_init, params_distance);

last_secrecy_rate = inf;
stall_count = 0;
stop_reason = "max_iter";
secrecy_rate_curve = [];
final_rho = nan;

for iter = 1:max_iter
    params_distance.theta_prev = v;
    params_distance.fixed_beamformer = f_current;
    obj_fun1 = @(x) fixed_v_fun(x, v, params_distance);

    initial_swarm = lhsdesign(num_particles, num_vars) .* (x_ub - x_lb) + x_lb;
    initial_swarm(1, :) = x_opt_distance;
    options = optimoptions('particleswarm', ...
        'SwarmSize', num_particles, ...
        'Display', 'off', ...
        'UseParallel', params_distance.pso_use_parallel, ...
        'FunctionTolerance', 1e-4, ...
        'MaxStallIterations', 20, ...
        'InitialSwarmMatrix', initial_swarm, ...
        'InertiaRange', [0.1, 0.9], ...
        'SocialAdjustmentWeight', 2, ...
        'SelfAdjustmentWeight', 1.6, ...
        'MaxIterations', 200);
    x_opt_distance = particleswarm(obj_fun1, num_vars, x_lb, x_ub, options);
    p_R = x_opt_distance(1:2).';
    psi = x_opt_distance(3:5).';

    params_distance.phase_surrogate_mode = true;
    [v, ~, ~] = rieman_grad_fun(p_R, psi, x_opt_distance, v, params_distance);
    if isfield(params_distance, 'phase_surrogate_mode')
        params_distance = rmfield(params_distance, 'phase_surrogate_mode');
    end

    [final_rho,~, hc_tmp, ht_tmp, ~, ~, ~, ~, ~, ~, ~] = bulid_H(p_R, psi, v, params_distance);
    f_tmp = solve_beamforming_vector(hc_tmp, ht_tmp, params_distance);
    params_distance.fixed_beamformer = f_tmp;
    [~, secrecy_rate_iter, ~, ~, ~] = secrecy_rate_objective(hc_tmp, ht_tmp, params_distance);
    secrecy_rate_curve = [secrecy_rate_curve; secrecy_rate_iter]; %#ok<AGROW>

    f_current = f_tmp;
    if isfield(params_distance, 'fixed_beamformer')
        params_distance = rmfield(params_distance, 'fixed_beamformer');
    end

    if isfinite(last_secrecy_rate)
        rel_progress = abs(secrecy_rate_iter - last_secrecy_rate) / max(abs(last_secrecy_rate), 1e-12);
    else
        rel_progress = inf;
    end
    if iter >= params_distance.min_ao_iter && rel_progress < params_distance.ao_rel_progress_tol
        stall_count = stall_count + 1;
    else
        stall_count = 0;
    end
    if stall_count >= params_distance.ao_stall_window
        stop_reason = "stall_window";
        break;
    end
    last_secrecy_rate = secrecy_rate_iter;
end

[is_converged, final_abs_gap, final_rel_gap] = certify_convergence(secrecy_rate_curve, params_distance);
if ~is_converged && stop_reason == "stall_window"
    stop_reason = "failed_final_certificate";
end

if isempty(secrecy_rate_curve)
    final_secrecy_rate = nan;
else
    final_secrecy_rate = secrecy_rate_curve(end);
end
end

function [is_converged, final_abs_gap, final_rel_gap] = certify_convergence(rate_curve, params_check)
%CERTIFY_CONVERGENCE 用最后一个窗口的真实保密率变化认证小实验是否已经稳定。
if isempty(rate_curve) || length(rate_curve) < 2
    is_converged = false;
    final_abs_gap = inf;
    final_rel_gap = inf;
    return;
end

window_len = min(params_check.convergence_check_window, length(rate_curve) - 1);
recent_rates = rate_curve(end-window_len:end);
recent_abs_gaps = abs(diff(recent_rates));
recent_rel_gaps = recent_abs_gaps ./ max(abs(recent_rates(1:end-1)), 1e-12);
final_abs_gap = max(recent_abs_gaps);
final_rel_gap = max(recent_rel_gaps);

% 收敛认证采用常见的“相对变化或绝对变化二选一”标准：
% 速率本身较大时看相对变化，速率接近0时看绝对变化，避免把已进入平台期的曲线误判为未收敛。
is_converged = length(rate_curve) >= params_check.min_ao_iter && ...
    (final_rel_gap < params_check.ao_rel_progress_tol || ...
    final_abs_gap < params_check.convergence_abs_tol);
end
