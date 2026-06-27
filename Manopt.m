clear; 
close all;
% 关闭有限差分警告（已提供精确解析海森矩阵）
warning('off', 'manopt:getHessian:approx'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================== 1. 参数设置与物理环境初始化 ==================
Flag.save = 0;
nx = 8; % 测试 4x4 网格，验证是否能跑出正确的 V 型或拱形
ny = 8;

%% Material constants & 物理矩阵装配
V0 = 0.1;
x_min = 1e-8;
scale0 = 1;
eeee = 20000*scale0;
rho = 7.86*10^(-4)*scale0;
dmm0 = 1*scale0;

% 调用外部组装函数
[dll, matH, coord_x, ir, irr, ird, Idx] = member(nx, ny, Inf);

nk = size(coord_x, 1);
nd = size(matH, 1); % 自由度 n
nm = size(matH, 2); % 杆件数 m

dmm = zeros(nd, 1);
dmm((nx-1):(nx), 1) = dmm0 * ones(2, 1); % 底部中心加载红点
ns_M = sparse(diag(dmm));

sqrtK = matH * sparse(diag(sqrt(eeee ./ dll)));
dm = elem_matrix(coord_x, irr, dll, rho, nm);
strm_M = mk_matrix(dm, ir, nm, nd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================== 2. 基于 Manopt 的黎曼流形 BM-ALM ==================
fprintf(' =========================================================== \n');
fprintf(' ==  Manopt BM-ALM (物理稳定性增强版) >>> \n');
fprintf('     #members = %g ;  #DOF = %g\n', nm, nd);

% 二分法区间设置
lambdaL = -1000;
lambdaU = 0;
maxiter_bisect = 25; 
tol_bisect = 1e-3;   

t_total = 0;
iter = 0;

% --- Manopt 流形环境初始化 ---
% 【修复1】针对小网格提升保底秩 (从 4 提升到 8)，确保二阶算法有足够维度逃逸局部最优
rank_r = min(18, max(8, round(nd / 4))); 
fprintf('     [参数自适应] 截断秩 rank_r 设定为: %d \n', rank_r);

manifold = symfixedrankYYfactory(nd, rank_r);
problem.M = manifold;
options.verbosity = 0; 
options.maxiter = 5;      
options.maxtime = 0.5;    
options.minstepsize = 1e-10; 

% 【修复2】大幅增强初始随机扰动 (从 0.1 提升到 0.5)，防止卡在绝对对称的中间柱子陷阱
x = (ones(nm, 1) + 0.5 * rand(nm, 1)); 
x = x * (V0 / (dll' * x)); 

V_opt = manifold.rand(); 
x_best = x; 

while abs(lambdaU - lambdaL) > tol_bisect && iter < maxiter_bisect
    tic
    lambda = (lambdaL + lambdaU) / 2;
    iter = iter + 1;
    
    U = zeros(nd, nd);
    rho_pen = 1.0;
    old_residual = Inf; 
    
    % ================= ALM 交替优化内层循环 =================
    for alm_iter = 1:200
        % ---------- Step 1: 求解 V (流形更新) ----------
        mat_K = sqrtK * sparse(diag(x)) * sqrtK';
        mat_M = reshape(strm_M * x, nd, nd);
        S_mat = mat_K + lambda * (mat_M + ns_M);
        S_mat = (S_mat + S_mat') / 2; 
        
        A_tgt = S_mat + U / rho_pen;
        A_tgt = (A_tgt + A_tgt') / 2;
        
        problem.cost  = @(V) 0.5 * norm(V*V' - A_tgt, 'fro')^2;
        problem.egrad = @(V) 2 * (V*V' - A_tgt) * V;
        problem.ehess = @(V, dV) 2 * ((dV*V' + V*dV')*V + (V*V' - A_tgt)*dV);
        
        [V_opt, ~, ~, ~] = trustregions(problem, V_opt, options);
        
        % ---------- Step 2: 求解 x (物理拓扑演化) ----------
        A_tgt_x = V_opt * V_opt' - U / rho_pen;
        for inner = 1:15 
            mat_K = sqrtK * sparse(diag(x)) * sqrtK';
            mat_M = reshape(strm_M * x, nd, nd);
            S_mat = mat_K + lambda * (mat_M + ns_M);
            S_mat = (S_mat + S_mat') / 2;
            
            W = S_mat - A_tgt_x; 
            grad_K = sum((W * sqrtK) .* sqrtK, 1)';
            grad_M = strm_M' * W(:);
            grad_x = grad_K + lambda * grad_M;
            
            % 【修复3】大幅降低单步最大变化率 (从 0.05 降至 0.01)
            % 这样杆件会慢慢向支座“伸展”，而不是瞬间堆积在中间形成孤立柱子
            max_change = 0.01 * V0; 
            grad_norm = norm(grad_x) + 1e-12;
            step_size = min(0.1, max_change / grad_norm); 
            x_new = x - step_size * grad_x;
            
            % 投影至体积约束面
            fun = @(mult) dll'*max(x_new - mult*dll, x_min) - V0;
            lb = min(x_new - V0)/max(dll) - 1e-6;
            ub = max(x_new)/min(dll) + 1e-6;
            mult = fzero(fun, [lb, ub]);
            x = max(x_new - mult*dll, x_min);
        end
        
        % ---------- Step 3: 更新惩罚项 ----------
        mat_K = sqrtK * sparse(diag(x)) * sqrtK';
        mat_M = reshape(strm_M * x, nd, nd);
        S_mat = mat_K + lambda * (mat_M + ns_M);
        S_mat = (S_mat + S_mat') / 2;
        
        residual = norm(S_mat - V_opt*V_opt', 'fro') / (norm(S_mat, 'fro') + 1e-6);
        
        if residual < 5e-3 
            break;
        end
        
        U = U + rho_pen * (S_mat - V_opt*V_opt');
        
        % 动态 rho 更新策略，保持惩罚参数平滑增长
        if residual > 0.25 * old_residual
            rho_pen = min(rho_pen * 1.5, 1e4); 
        end
        old_residual = residual;
    end
    % =======================================================
    
    if residual < 5e-3
        lambdaU = lambda; 
        flag = 0;
        x_best = x;       
        
        figure(1); 
        [~] = draw_cs_freq(coord_x, irr, x_best, 1, nx, ny);
        drawnow; 
    else
        lambdaL = lambda; 
        flag = 1;
    end
    
    t_iter = toc;
    t_total = t_total + t_iter;
    fprintf('Iter.:%4i Obj.:%15.10f Feas:%1i ALM_Iters:%3i Time:%8.2f \n', ...
        iter, lambda, flag, alm_iter, t_total); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================== 3. Output 结果输出 ==================
fprintf(' =========================================================== \n');
fprintf(' 最终特征值 = %15.10f \n', lambdaU);
fprintf(' =========================================================== \n');

figure(1);
[~] = draw_cs_freq(coord_x, irr, x_best, 1, nx, ny);
title(sprintf('Final Topology (Grid: %dx%d, Freq: %.4f)', nx, ny, -lambdaU));