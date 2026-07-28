% function [v_opt, f_value, info] = rieman_grad_fun(p_R, psi, x_opt, v, params)
% import manopt.framework.*   % 导入 Manopt 框架
% import manopt.solvers.*     % 导入求解器
% % 定义球面流形
% n = params.N;  %变量的维度
% M = complexcirclefactory(n); %定义复球面流形，每个分量的模值为1
% 
% % 目标函数
% obj_fun2 = @(v)fixed_pR_psi_fun(x_opt, v, params);
% euclid_grad = @(v) euclid_grad_fun(p_R, psi, v, params);  % 欧式梯度
% 
% %定义 problem 结构体 (确保包含 cost 和 grad) 
% problem.M = M;
% problem.cost = @(v) obj_fun2(v);
% problem.egrad = @(v) euclid_grad (v);  
% problem.grad = @(v) M.egrad2rgrad(v, problem.egrad(v));  % 计算黎曼梯度
% 
% % 使用梯度下降法更新v，并投影回球面(使用 Manopt 提供的优化器)
% options.verbosity = 1;  %表示设定输出级别为 2，即每次迭代都会显示优化过程的基本信息。
% options.maxiter = 2000; % 最大迭代次数
% options.tolgradnorm = 1e-10; %梯度范数容差
% options.minstepsize = 1e-10;
% [v_opt, f_value, info] = steepestdescent(problem, v, options);  % 运行梯度下降优化
% % 信赖域方法：trustregions
% % 最速下降法：steepestdescent
% end

function [v_opt, f_value, info] = rieman_grad_fun(p_R, psi, x_opt, v, params)
import manopt.framework.*   % 导入 Manopt 框架
import manopt.solvers.*     % 导入求解器
%RIEMAN_GRAD_FUN Phase update with an inner MM loop.
% Input/output v is the physical IRS diagonal phi=diag(Phi). The paper phase
% variable is theta=conj(phi), because Phi=diag(theta^*). Therefore Manopt
% optimizes theta internally and converts back to phi at the end.
% For fixed g and beamformer, repeat theta_k -> beta_k/lambda_k -> theta_{k+1}.

n = params.N;
M = complexcirclefactory(n);

if isfield(params, 'phase_mm_maxiter')
    phase_mm_maxiter = params.phase_mm_maxiter;
else
    phase_mm_maxiter = 10;
end
if isfield(params, 'phase_mm_tol')
    phase_mm_tol = params.phase_mm_tol;
else
    phase_mm_tol = 1e-4;
end

% Manopt options for each surrogate minimization.
if isfield(params, 'manopt_verbosity')
    options.verbosity = params.manopt_verbosity;
else
    options.verbosity = 1;
end
if isfield(params, 'manopt_maxiter')
    options.maxiter = params.manopt_maxiter;
else
    options.maxiter = 2000;
end
if isfield(params, 'manopt_tolgradnorm')
    options.tolgradnorm = params.manopt_tolgradnorm;
else
    options.tolgradnorm = 1e-10;
end
if isfield(params, 'manopt_minstepsize')
    options.minstepsize = params.manopt_minstepsize;
else
    options.minstepsize = 1e-10;
end

theta_current = conj(v);
mm_step_norm = nan(phase_mm_maxiter, 1);
mm_cost = nan(phase_mm_maxiter, 1);
mm_inner_info = cell(phase_mm_maxiter, 1);
f_value = nan;

for mm_iter = 1:phase_mm_maxiter
    % Rebuild beta/lambda around the current MM point theta_prev.
    local_params = params;
    local_params.theta_prev = theta_current;

    % fixed_pR_psi_fun/bulid_H use the physical phase phi, while Manopt uses theta.
    obj_fun2 = @(theta) fixed_pR_psi_fun(x_opt, conj(theta), local_params);
    euclid_grad = @(theta) euclid_grad_fun(p_R, psi, theta, local_params);

    problem.M = M;
    problem.cost = @(theta) obj_fun2(theta);
    problem.egrad = @(theta) euclid_grad(theta);
    problem.grad = @(theta) M.egrad2rgrad(theta, problem.egrad(theta));

    [theta_next, f_value, inner_info] = steepestdescent(problem, theta_current, options);
    mm_inner_info{mm_iter} = inner_info;
    mm_cost(mm_iter) = f_value;
    mm_step_norm(mm_iter) = norm(theta_next - theta_current, 2) / sqrt(n);

    theta_current = theta_next;
    if mm_step_norm(mm_iter) < phase_mm_tol
        break;
    end
end

v_opt = conj(theta_current);
info = struct;
info.mm_iterations = mm_iter;
info.mm_step_norm = mm_step_norm(1:mm_iter);
info.mm_cost = mm_cost(1:mm_iter);
info.inner_info = mm_inner_info(1:mm_iter);
info.stop_reason = 'max_phase_mm_iter';
if mm_step_norm(mm_iter) < phase_mm_tol
    info.stop_reason = 'phase_mm_step_tol';
end
end
