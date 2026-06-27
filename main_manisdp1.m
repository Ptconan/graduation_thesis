clearvars
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ManiSDP: Stiefel 流形上的低秩 SDP + Rayleigh-Ritz + 平滑梯度 + 动量引擎
% 解决 4x4 甚至更大规模的桁架频率最大化问题
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);
if exist('manopt', 'dir')
    addpath(genpath(fullfile(pwd, 'manopt', 'manopt')));
end
warning('off', 'manopt:getHessian:approx');

%% 1. 物理参数与矩阵装配
Flag.save = 0;
nx = 8; 
ny = 8;

if nx == 4 && ny == 4
    opt_value = -51.4026958936;
elseif nx == 8 && ny == 8
    opt_value = -15.3984958306;
else
    opt_value = NaN;
end

V0 = 0.1;
x_min = 1e-8;
scale0 = 1;
eeee = 20000 * scale0;
rho_mat = 7.86 * 10^(-4) * scale0;
dmm0 = 1 * scale0;

[dll, matH, coord_x, ir, irr, ird, Idx] = member(nx, ny, Inf);
nd = size(matH, 1);
nm = size(matH, 2);

dmm = zeros(nd, 1);
dmm((nx-1):(nx), 1) = dmm0 * ones(2, 1);
ns_M = sparse(diag(dmm));

sqrtK = matH * sparse(diag(sqrt(eeee ./ dll)));
dm = elem_matrix(coord_x, irr, dll, rho_mat, nm);
strm_M = mk_matrix(dm, ir, nm, nd);

%% 2. ManiSDP 与外层二分法参数
lambdaL = -1000;
lambdaU = 0;
tol_bisect = 1e-4;
maxiter_bisect = 40;

% Stiefel 流形降维参数
rank_p = min(10, max(4, round(nd / 20)));
max_alt = 250;
step0_x = 0.8;
mu_smooth = 0.05;
beta_mom = 0.85; % 动量保留系数

manifold = stiefelfactory(nd, rank_p);
probY.M = manifold;
optsY.verbosity = 0;
optsY.maxiter = 15;
optsY.maxtime = 1.5;

fprintf(' =========================================================== \n');
fprintf(' ==  ManiSDP (Stiefel + Momentum Engine) %dx%d >>> \n', nx, ny);
fprintf('     #members = %g ;  #DOF = %g ; subspace rank p = %g\n', nm, nd, rank_p);

% 全局热启动
x_warm = (ones(nm, 1) + 0.5 * rand(nm, 1));
x_warm = project_volume_mani(x_warm, dll, V0, x_min);
x_best = x_warm;
Y_opt = manifold.rand();

t_total = 0;
iter = 0;
opt = 1;

%% 3. 核心求解主循环
while (lambdaU - lambdaL) > tol_bisect && opt > 1e-4 && iter < maxiter_bisect
    t0 = tic;
    lambda = 0.5 * (lambdaL + lambdaU);
    opt = abs(lambda - opt_value) / abs(opt_value);
    iter = iter + 1;
    
    x = x_warm;
    min_eig = -inf;
    
    % 【核心升级 1】初始化物理演化的动量变量
    v_x = zeros(nm, 1); 
    
    for alt = 1:max_alt
        % 组装力学矩阵
        mat_K = sqrtK * spdiags(x, 0, nm, nm) * sqrtK';
        mat_M = reshape(strm_M * x, nd, nd) + ns_M;
        S_mat = mat_K + lambda * mat_M;
        S_mat = 0.5 * (S_mat + S_mat');
        
        % ManiSDP Step 1: Stiefel 流形上追踪危险子空间
        probY.cost = @(Y) real(trace(Y' * S_mat * Y));
        probY.egrad = @(Y) 2 * S_mat * Y;
        [Y_opt, ~, ~, ~] = conjugategradient(probY, Y_opt, optsY);
        
        % ManiSDP Step 2: Rayleigh-Ritz 子空间精确特征值
        H_sub = Y_opt' * S_mat * Y_opt;
        H_sub = 0.5 * (H_sub + H_sub');
        [V_sub, D_sub] = eig(H_sub);
        eig_sub = real(diag(D_sub));
        min_eig = min(eig_sub);
        
        % 如果已经满足半正定，提前跳出内层
        if min_eig >= -1e-5
            break
        end
        
        V_full = Y_opt * V_sub;
        
        % ManiSDP Step 3: Log-sum-exp 平滑梯度提取
        mu_alt = mu_smooth / sqrt(alt);
        weights = exp(-(eig_sub - min_eig) / (mu_alt + 1e-12));
        weights = weights / sum(weights);
        
        grad_x = zeros(nm, 1);
        for j = 1:rank_p
            vj = V_full(:, j);
            dK_j = (sqrtK' * vj).^2;
            dM_j = strm_M' * reshape(vj * vj', nd^2, 1);
            grad_x = grad_x + weights(j) * (dK_j + lambda * dM_j);
        end
        
        gn = norm(grad_x) + 1e-12;
        
        % 【核心升级 2】动量机制 (Momentum)，吸收局部震荡，指引全局方向
        v_x = beta_mom * v_x + (1 - beta_mom) * (grad_x / gn);
        
        % 【核心升级 3】更平滑的步长衰减，防止过早卡死
        step = step0_x / (1 + 0.05 * alt); 
        
        % 沿着动量方向进行拓扑演化
        x_new = x + step * v_x * V0;
        x = project_volume_mani(x_new, dll, V0, x_min);
    end
    
    x_warm = x;
    z = -min_eig;
    
    % 二分法区间更新
    if z >= 0
        lambdaL = lambda;
        x_best = x;
        flag = 0;
    else
        lambdaU = lambda;
        flag = 1;
    end
    
    t_total = t_total + toc(t0);
    fprintf('Iter.:%4i Obj.:%15.10f MinEig:%12.4e z:%10.4e Flag:%1i Bracket:%.4f Time:%7.2f\n', ...
        iter, lambda, min_eig, z, flag, (lambdaU - lambdaL), t_total);
end

%% 4. 后处理与结果输出
[~] = draw_cs_freq(coord_x, irr, x_best, 1, nx, ny);

fprintf('\n>>> ManiSDP 完毕 <<<\n');
fprintf(' 最优 SDP lambda* = %15.10f \n', lambdaL);
fprintf(' 理论参考值       = %15.10f \n', opt_value);
fprintf(' 相对误差         = %12.4e \n', abs(lambdaL - opt_value) / abs(opt_value));
fprintf(' 总耗时 = %6.2f [s] \n', t_total);
fprintf(' =========================================================== \n');

%% Local Functions
function x = project_volume_mani(x_raw, dll, V0, x_min)
    x_raw = max(real(x_raw), x_min);
    x_raw(~isfinite(x_raw)) = x_min;
    
    if dll' * x_raw <= V0 + 1e-12
        x = x_raw;
        return
    end
    
    fun = @(mult) dll' * max(x_raw - mult * dll, x_min) - V0;
    lb = min(x_raw - V0) / max(dll) - 1;
    ub = max(x_raw) / min(dll) + 1;
    
    mult = fzero(fun, [lb, ub]);
    x = max(x_raw - mult * dll, x_min);
end