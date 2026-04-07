% 自定义迭代输出函数
function stop = pso_output_function(optimValues, state)
    stop = false; % 默认不停止优化
    
    % 仅在 "iter" 状态下打印信息
    if strcmp(state, 'iter')
        fprintf('迭代 %d: 最优函数值 = %.6f\n', ...
            optimValues.iteration, optimValues.bestfval);
    end
end