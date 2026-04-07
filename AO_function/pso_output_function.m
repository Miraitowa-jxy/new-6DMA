% % 自定义迭代输出函数
% function stop = pso_output_function(optimValues, state)
%     stop = false; % 默认不停止优化
% 
%     % 仅在 "iter" 状态下打印信息
%     if strcmp(state, 'iter')
%         fprintf('迭代 %d: 最优函数值 = %.6f\n', ...
%             optimValues.iteration, optimValues.bestfval);
%     end
% end
% 

function stop = pso_output_function(optimValues, state)
    % 声明全局变量
    global bestval_history
    stop = false; % 默认不停止优化

    if strcmp(state, 'init')
        bestval_history = [];  % 初始化变量
    elseif strcmp(state, 'iter')
        % 记录每一次的最优函数值
        bestval_history = [bestval_history; optimValues.bestfval];

        % 也可以输出当前迭代信息
        fprintf('迭代 %d: 最优函数值 = %.6f\n', ...
            optimValues.iteration, optimValues.bestfval);
    end
end
