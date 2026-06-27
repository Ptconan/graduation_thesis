% =========================================================================
% 论文算法：基于 Burer-Monteiro 低秩分解与增广拉格朗日法 (BM-ALM) 的桁架拓扑优化
% 功能：求解大规模桁架结构在体积约束下的特征值最大化等价问题
% =========================================================================

clearvars; close all; clc;
tic;

%% 1. 物理参数与环境初始化 (完全对齐 main_bisection.m)
nx = 8; ny = 8;           % 网格规模 (可修改为 8x8 测试大规模优势)
V0 = 0.1;                 % 体积分数上限
x_min = 1e-8;             % 杆件截面积下界
scale0 = 1;
eeee = 20000 * scale0;    % 杨氏模量 E
rho_material = 7.86e-4;   % 材料密度
dmm0 = 1 * scale0;        % 非结构质量强度

% 目标特征值 (此处填入内点法算出的真值，或设为一个较高的目标进行可行性测试)
% 对于 4x4, Inf 情况，理论最优 lambda 约为 51.4026
lambda_target = 15.3984958306; 

% 调用基结构生成函数 (确保 member.m 在同一路径下)
[dll, matH, coord_x, ir, irr, ird, Idx] = member(nx, ny, Inf);
nd = size(matH, 1);       % 系统自由度 (N)
nm = size(matH, 2);       % 杆件总数 (m)

% 构建非结构质量矩阵 M0
dmm = zeros(nd, 1);
dmm((nx-1):nx, 1) = dmm0 * ones(2, 1);
ns_M = sparse(diag(dmm));

%% 2. BM-ALM 算法特定参数
rank_r = 10;               % 低秩矩阵的秩 (r << N)
rho_alm = 1.0;            % 初始增广拉格朗日惩罚参数
max_iter = 500;           % 最大迭代次数
tol = 1e-4;               % 收敛容差

fprintf('-----------------------------------------------------------\n');
fprintf('Initializing Low-Rank BM-ALM Framework\n');
fprintf('Grid: %dx%d | Variables: %d | Rank: %d\n', nx, ny, nm, rank_r);
fprintf('-----------------------------------------------------------\n');

%% 3. 变量初始化
x = (V0 / sum(dll)) * ones(nm, 1); % 初始杆件均匀分布
V = randn(nd, rank_r) * 0.01;      % 初始低秩矩阵 V (维度 N x r)
Y = zeros(nd, nd);                 % 拉格朗日乘子矩阵 (对偶变量)

% 存储收敛数据用于绘图
res_history = zeros(max_iter, 1);

%% 4. 优化主循环
for iter = 1:max_iter
    
    % --- 步骤 A: 组装当前物理刚度矩阵 K(x) ---
    % 注意：桁架单位刚度矩阵 K_i = (E/L_i) * h_i * h_i'
    K_x = eeee * matH * diag(x ./ dll) * matH'; 
    M_x = ns_M; % 本算例主要考虑非结构质量
    
    % A(x) 为优化问题的核心矩阵项
    A_x = K_x - lambda_target * M_x; 
    
    % --- 步骤 B: 更新低秩矩阵 V (Burer-Monteiro 分解步) ---
    % 目标：让 V*V' 逼近 (A_x + Y/rho)
    Target_Mat = A_x + Y / rho_alm;
    Target_Mat = (Target_Mat + Target_Mat') / 2; % 保持对称
    
    % 特征值分解获取最优低秩逼近
    [EigVec, EigVal] = eigs(Target_Mat, rank_r, 'largestreal'); 
    EigVal_diag = diag(EigVal);
    EigVal_diag(EigVal_diag < 0) = 0; % 投影到半正定锥
    V = EigVec * diag(sqrt(EigVal_diag));
    
    % --- 步骤 C: 更新截面积 x (归一化梯度投影步) ---
    % 计算增广拉格朗日函数关于 x 的梯度
    % 这里的梯度计算利用了 A_x 关于 x 的线性关系
    Error_Mat = A_x - V * V' + Y / rho_alm;
    
    grad_x = zeros(nm, 1);
    for i = 1:nm
        h_i = matH(:, i);
        % 复合导数推导：d(L)/dx_i
        grad_x(i) = (eeee/dll(i)) * (h_i' * Error_Mat * h_i); 
    end
    
    % 梯度归一化：消除物理量纲（如20000）带来的数值不稳定
    grad_x = grad_x / (norm(grad_x) + 1e-12);
    
    % 执行投影梯度下降
    step_size = (V0 / 5) / sqrt(iter); % 动态衰减步长
    x_trial = x - step_size * grad_x;
    
    % 鲁棒体积投影 (Bisection search on mu)
    mu_L = -1e8; mu_U = 1e8;
    for b_iter = 1:50
        mu_mid = (mu_L + mu_U) / 2;
        x_proj = max(x_trial - mu_mid * dll, x_min);
        if sum(x_proj .* dll) > V0
            mu_L = mu_mid;
        else
            mu_U = mu_mid;
        end
    end
    x = x_proj;
    
    % --- 步骤 D: 更新拉格朗日乘子 Y 与 惩罚参数 rho ---
    % 重新组装 K(x) 以更新残差
    K_x_updated = eeee * matH * diag(x ./ dll) * matH';
    Residual = (K_x_updated - lambda_target * M_x) - V * V';
    
    % 乘子更新（ALM 核心逻辑）
    Y = Y + rho_alm * Residual;
    
    % 动态提升惩罚系数以加速收敛
    rho_alm = min(rho_alm * 1.05, 1000); 
    
    % --- 数据记录与日志 ---
    res_norm = norm(Residual, 'fro');
    res_history(iter) = res_norm;
    
    if mod(iter, 10) == 0 || iter == 1
        fprintf('Iter: %3d | Residual: %.4e | rho: %.2f\n', iter, res_norm, rho_alm);
    end
    
    % 收敛判定
    if res_norm < tol && iter > 20
        fprintf('Algorithm converged at iteration %d!\n', iter);
        break;
    end
end

total_time = toc;
fprintf('-----------------------------------------------------------\n');
fprintf('Optimization Finished!\n');
fprintf('Total CPU Time: %.4f seconds\n', total_time);
fprintf('Final Residual: %.4e\n', res_norm);
fprintf('-----------------------------------------------------------\n');

%% 5. 结果可视化
% 绘制残差收敛曲线
figure('Color', 'w');
semilogy(res_history(1:iter), 'LineWidth', 2);
grid on; xlabel('Iteration'); ylabel('Constraint Residual (Log Scale)');
title('Convergence of Low-Rank BM-ALM');

% 绘制最终拓扑结构图 (调用 draw_cs_freq)
if exist('draw_cs_freq', 'file')
    draw_cs_freq(coord_x, irr, x, 10, nx, ny);

else
    fprintf('Warning: draw_cs_freq.m not found. Visual map skipped.\n');
end