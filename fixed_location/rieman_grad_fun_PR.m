function [v_opt, f_value, info] = rieman_grad_fun_PR(p_R, psi, x_opt, v, params)
import manopt.framework.*   % 导入 Manopt 框架
import manopt.solvers.*     % 导入求解器
% 定义球面流形
n = params.N;  %变量的维度
M = complexcirclefactory(n); %定义复球面流形，每个分量的模值为1

% 目标函数
obj_fun2 = @(v)fixed_pR_psi_fun_PR(x_opt, v, params);
euclid_grad = @(v) euclid_grad_fun_PR(p_R, psi, v, params);  % 欧式梯度

%定义 problem 结构体 (确保包含 cost 和 grad) 
problem.M = M;
problem.cost = @(v) obj_fun2(v);
problem.egrad = @(v) euclid_grad (v);  
problem.grad = @(v) M.egrad2rgrad(v, problem.egrad(v));  % 计算黎曼梯度

% % 使用梯度下降法更新v，并投影回球面(使用 Manopt 提供的优化器)
% options.verbosity = 1;  %表示设定输出级别为 2，即每次迭代都会显示优化过程的基本信息。
% options.maxiter = 2000; % 最大迭代次数
% options.tolgradnorm = 1e-6; %梯度范数容差
% options.minstepsize = 1e-6;

%不优化姿态的话必须使用这个才能保证收敛,控制V变量收敛，保证曲线平滑性，任何时候不平滑都可调
options.verbosity = 1;  
options.maxiter = 2000; 
options.tolgradnorm = 1e-10; 
options.minstepsize = 1e-10;

[v_opt, f_value, info] = steepestdescent(problem, v, options);  % 运行梯度下降优化
% 信赖域方法：trustregions
% 最速下降法：steepestdescent
end