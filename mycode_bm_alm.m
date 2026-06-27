%%%%% 结构拓扑优化 BM-ALM + minFunc (满秩终极版)
%%%%% 适用于复现大规模结构优化的基准测试
clear; clc; % 全程序只允许在此处清理一次工作区

% --- 1. 参数初始化与数据加载 ---
nd = 46;           % 自由度
n = 200;           % 杆件截面积变量数
full_nd = 2*nd+1;  % 矩阵总维度 (93)
lambda = 5e6;

load('sqrtK.txt');
load('strmM.txt');
fl = zeros(nd, 1); fl(41) = 4;

% 组装常数项 A0
disp('正在组装系统特征矩阵...');
A0 = [1, -fl', zeros(1,nd); 
     -fl, zeros(nd,nd), zeros(nd,nd); 
     zeros(nd,1), zeros(nd,nd), zeros(nd,nd)];
 
% 高度向量化 A_col (d^2 x n)，单次组装彻底消灭后续 for 循环
A_col = zeros(full_nd * full_nd, n);
for i = 1:n
    Ai = [zeros(1, 2*nd+1); 
          zeros(nd, 1), sqrtK(:,i)*sqrtK(:,i)', zeros(nd, nd); 
          zeros(nd, nd+1), sqrtK(:,i)*sqrtK(:,i)' - lambda*reshape(strmM(:,i), nd, nd)];
    A_col(:, i) = Ai(:);
end

% --- 2. 算法初始化 (终极满秩设定) ---
p = full_nd;                          % 【核心修改】令 p = 93 满秩表达，彻底容纳结构受力的满秩矩阵！
V = randn(full_nd, p) * 0.01;         % 初始化正方形 V 矩阵 (93x93)
y = ones(n, 1);                       % 初始化 y (此时截面积 x = y.^2 = 1)
Y = zeros(full_nd, full_nd);          % 拉格朗日乘子矩阵
rho = 1.0;                            % 初始惩罚项 
tol = 1e-4;                           % 目标收敛残差容差

% --- 3. 配置 minFunc 求解器参数 ---
options_mf = [];
options_mf.method = 'lbfgs';      
options_mf.display = 'none';          % 关闭内层刷屏
options_mf.maxIter = 500;             % 内层给予充足迭代次数找准梯度方向
options_mf.DerivativeCheck = 'off';

disp('>> 满秩空间 BM-ALM + minFunc 引擎全速启动 <<');
tic;

% --- 4. ALM 外层主循环 ---
% 满秩由于搜索空间大，我们给予外层最多 60 次迭代机会
for outer_iter = 1:60
    % 拼接当前的无约束决策变量
    z0 = [y; V(:)];
    
    % 调用一阶极速引擎 minFunc 求解无约束子问题
    z_opt = minFunc(@(z) alm_loss_gradient(z, A0, A_col, Y, rho, n, full_nd, p), z0, options_mf);
    
    % 解耦优化结果
    y = z_opt(1:n);
    V = reshape(z_opt(n+1:end), full_nd, p);
    
    % 恢复拓扑截面积变量 x (隐式满足 x >= 0)
    x = y.^2;
    
    % 计算当前 LMI 约束的违例程度 (残差)
    Ax = A0 + reshape(A_col * x, full_nd, full_nd);
    residual = Ax - V * V';
    res_norm = norm(residual, 'fro');
    
    % 打印规范的科研收敛日志
    fprintf('外层 Iter %2d | 结构体积: %8.4f | 约束残差: %.2e | 惩罚项 rho: %.1f\n', ...
            outer_iter, sum(x), res_norm, rho);
        
    % 完美达到收敛精度，提前退出
    if res_norm < tol
        fprintf('\n>> 恭喜！算法完美达到预设精度(%.2e)，成功收敛！<<\n', tol);
        break;
    end
    
    % 乘子与惩罚项标准更新
    Y = Y + rho * residual;
    rho = min(rho * 1.5, 1e7);        % 适当加快一点惩罚放大速度逼迫残差收敛
end
toc;

% 最终的最优拓扑截面积保存在变量 x_opt 中
x_opt = x; 
fprintf('\n最终结构总体积: %.4f (参考基准 CVX: 922.929)\n', sum(x_opt));

% 尝试绘制结果 (如果当前目录下有你之前的绘图脚本)
try
    [~,~,coord_x,~,irr,~,~] = member(4,4,Inf);
    draw_cs_freq(coord_x, irr, x_opt, 1, 4, 4);
    disp('绘图成功。');
catch
    disp('未检测到绘图函数，已跳过绘图。');
end


% =========================================================================
%  子函数：计算无约束光滑目标函数及其解析梯度 (必须放在脚本最末尾)
% =========================================================================
function [f, g] = alm_loss_gradient(z, A0, A_col, Y, rho, n, d, p)
    % 拆解变量
    y = z(1:n);
    V = reshape(z(n+1:end), d, p);
    x = y.^2;
    
    % 快速计算总矩阵 A(x)
    Ax = A0 + reshape(A_col * x, d, d);
    R = Ax - V * V';
    R_vec = R(:);
    
    % 1. 增广拉格朗日函数值
    f = sum(x) + Y(:)' * R_vec + (rho/2) * (R_vec' * R_vec);
    
    % 2. 梯度计算 (包含边界约束消除的链式法则导数传递)
    M_grad = Y + rho * R;  
    gy = 2 .* y .* (1 + A_col' * M_grad(:)); % 对 y 导数的链式法则修正
    gV = -2 * M_grad * V;
    
    g = [gy; gV(:)];
end