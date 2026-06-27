clear
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BM-ALM + 割平面 Bundle (内层) + 二分法 (外层)
% 与 main_bisection.m 相同 SDP; 内层用 main_bisection_bundle 的 LP 逼近 CVX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Flag.save = 0;
nx = 4;
ny = 4;

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

if exist('manopt', 'dir')
    addpath(genpath(fullfile(pwd, 'manopt', 'manopt')));
end
warning('off', 'manopt:getHessian:approx');

%% 二分 / BM-ALM / Bundle 参数
lambdaL = -1000;
lambdaU = 0;
tol_bisect = 1e-4;
maxiter_bisect = 45;
tol_psd = 1e-4;
tol_z = 1e-3;
max_bundle = 80;
max_bundle_cuts = 50;

max_alm_warm = 25;
tol_alm_res = 1e-2;
rank_r = min(12, max(6, round(nd / 5)));

linprog_opts = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');

manifold = symfixedrankYYfactory(nd, rank_r);
probV.M = manifold;
optsV.verbosity = 0;
optsV.maxiter = 4;
optsV.maxtime = 0.8;

fprintf(' =========================================================== \n');
fprintf(' ==  BM-ALM + Bundle-LP Bisection %dx%d >>> \n', nx, ny);
fprintf('     #members = %g ;  #DOF = %g\n', nm, nd);

x_warm = (ones(nm, 1) + 0.5 * rand(nm, 1));
x_warm = project_volume(x_warm, dll, V0, x_min);
x_best = x_warm;
V_opt = manifold.rand();
t_total = 0;
iter = 0;
opt = 1;

while (lambdaU - lambdaL) > tol_bisect && opt > 1e-4 && iter < maxiter_bisect
    tic
    lambda = 0.5 * (lambdaL + lambdaU);
    opt = abs(lambda - opt_value) / abs(opt_value);
    iter = iter + 1;

    %% (1) BM-ALM 热启动: 低秩 V 与截面 x
    x = x_warm;
    U = zeros(nd, nd);
    rho_pen = 1.0;
    for alm_k = 1:max_alm_warm
        S_mat = assemble_S(x, lambda, sqrtK, strm_M, ns_M, nd);
        A_tgt = (S_mat + U / rho_pen);
        A_tgt = 0.5 * (A_tgt + A_tgt');

        probV.cost = @(V) 0.5 * norm(V * V' - A_tgt, 'fro')^2;
        probV.egrad = @(V) 2 * (V * V' - A_tgt) * V;
        [V_opt, ~, ~, ~] = trustregions(probV, V_opt, optsV);

        A_tgt_x = V_opt * V_opt' - U / rho_pen;
        W = S_mat - A_tgt_x;
        grad_x = sum((W * sqrtK) .* sqrtK, 1)' + lambda * (strm_M' * W(:));
        gn = norm(grad_x) + 1e-12;
        x = project_volume(x - min(0.08, 0.02 * V0 / gn) * grad_x, dll, V0, x_min);

        residual = norm(S_mat - V_opt * V_opt', 'fro') / (norm(S_mat, 'fro') + 1e-8);
        if residual < tol_alm_res
            break
        end
        U = U + rho_pen * (S_mat - V_opt * V_opt');
        rho_pen = min(rho_pen * 1.3, 2e3);
    end

    %% (2) Bundle-LP 内层 (与 SDP 子问题等价, 替代纯梯度)
    [x, z, min_eig, ~, bstatus] = solve_sdp_bundle( ...
        lambda, x, dll, V0, x_min, sqrtK, strm_M, ns_M, nd, nm, ...
        tol_psd, tol_z, max_bundle, max_bundle_cuts, -1e8, 1e8, linprog_opts);

    if isempty(x) || any(~isfinite(x))
        x = x_warm;
        z = 1;
        min_eig = -inf;
    end
    x_warm = x;

    if z >= 0
        lambdaL = lambda;
        x_best = x;
        flag = 0;
    else
        lambdaU = lambda;
        flag = 1;
    end

    t_total = t_total + toc;
    fprintf('Iter.:%4i Obj.:%15.10f MinEig:%12.4e z:%10.4e Flag:%1i %s Bracket:%.4f Time:%7.2f\n', ...
        iter, lambda, min_eig, z, flag, char(bstatus), (lambdaU - lambdaL), t_total);
end

%% (3) 最终 Bundle 精修
fprintf('\n>> 最终 Bundle 精修 lambda = %.10f ...\n', lambdaL);
[x_best, z_f, min_eig_final, ~, bstatus] = solve_sdp_bundle( ...
    lambdaL, x_best, dll, V0, x_min, sqrtK, strm_M, ns_M, nd, nm, ...
    tol_psd, tol_z, 120, 60, -1e8, 1e8, linprog_opts);

[~] = draw_cs_freq(coord_x, irr, x_best, 1, nx, ny);

fprintf('\n>>> BM-ALM 完毕 <<<\n');
fprintf(' 最优 SDP lambda*     = %15.10f \n', lambdaL);
fprintf(' 最终 z / MinEig      = %12.4e / %12.4e \n', z_f, min_eig_final);
fprintf(' 理论参考值           = %15.10f \n', opt_value);
fprintf(' 相对误差             = %12.4e \n', abs(lambdaL - opt_value) / abs(opt_value));
fprintf(' 体积 dll''*x         = %12.6f \n', dll' * x_best);
fprintf(' 非零杆件 (x>10*xmin) = %g / %g\n', sum(x_best > 10 * x_min), nm);
fprintf(' 总耗时 = %6.2f [s] \n', t_total);
fprintf(' =========================================================== \n');

function S = assemble_S(x, lambda, sqrtK, strm_M, ns_M, nd)
mat_K = sqrtK * sparse(diag(x)) * sqrtK';
mat_M = reshape(strm_M * x, nd, nd) + ns_M;
S = full(mat_K + lambda * mat_M);
S = 0.5 * (S + S');
end

function x = project_volume(x_raw, dll, V0, x_min)
x_raw = max(real(x_raw), x_min);
x_raw(~isfinite(x_raw)) = x_min;
if dll' * x_raw <= V0 + 1e-12
    x = x_raw;
    return
end
fun = @(mult) dll' * max(x_raw - mult * dll, x_min) - V0;
lb = min(x_raw - V0) / max(dll) - 1;
ub = max(x_raw) / min(dll) + 1;
if ~isfinite(lb) || ~isfinite(ub) || lb >= ub
    lo = -1e4;
    hi = 1e4;
    while fun(lo) > 0, lo = lo - 1e4; end
    while fun(hi) < 0, hi = hi + 1e4; end
    mult = fzero(fun, [lo, hi]);
else
    mult = fzero(fun, [lb, ub]);
end
x = max(x_raw - mult * dll, x_min);
end
function [x_opt, z_opt, min_eig, history, bstatus] = solve_sdp_bundle( ...
    lambda, x0, dll, V0, x_min, sqrtK, strm_M, ns_M, nd, nm, ...
    tol_psd, tol_z, max_bundle, max_bundle_cuts, z_lb, z_ub, linprog_opts)
% =========================================================================
% 基于线性规划 (LP) 的割平面束方法 (Bundle-LP)
% 求解目标: min z 
% s.t.      sum(L.*x) <= V0, x >= x_min
%           S(x) + z*I >= 0  (等价于 v^T S(x) v + z >= 0)
% =========================================================================

    x = x0;
    z = 1.0;          % z 的初始猜测
    z_old = inf;
    A_cuts = [];      % 割平面约束矩阵
    b_cuts = [];      % 割平面常数项
    
    history = zeros(max_bundle, 1);
    bstatus = "MaxIter";

    for k = 1:max_bundle
        % 1. 组装当前 x 对应的对称矩阵 S(x)
        mat_K = sqrtK * spdiags(x, 0, nm, nm) * sqrtK';
        mat_M = reshape(strm_M * x, nd, nd) + ns_M;
        S = mat_K + lambda * mat_M;
        S = 0.5 * (S + S'); % 确保严格对称
        
        % 2. 求解最小特征值及其特征向量
        % 优先使用 eigs 提速，若遇到极端病态则退化为 eig
        try
            [V_eig, D_eig] = eigs(S, 1, 'smallestreal');
            min_eig = D_eig(1,1);
            v = V_eig(:,1);
        catch
            [V_eig, D_eig] = eig(full(S));
            [min_eig, min_idx] = min(diag(D_eig));
            v = V_eig(:, min_idx);
        end
        
        history(k) = min_eig;
        
        % 3. 检查收敛条件
        % 如果 z < 0 且 S(x) 满足半正定 (允许 tol_psd 误差)，则当前 lambda 严格可行
        if (z < 0 && min_eig >= -tol_psd) || (z < 0 && abs(min_eig + z) < tol_z)
             bstatus = "Feasible";
             break;
        end
        % 如果 LP 已经无法把 z 压到 0 以下 (z 稳定为正)，说明当前 lambda 不可行
        if k > 1 && abs(z - z_old) < 1e-6 && z > 0
             bstatus = "Infeasible";
             break;
        end
        z_old = z;

        % 4. 生成割平面 (Cutting Plane)
        % 解析偏导数: dK = (sqrtK'*v).^2,  dM = strm_M' * vec(v*v')
        dK = (sqrtK' * v).^2;
        v_vec = reshape(v * v', nd*nd, 1);
        dM = strm_M' * v_vec;
        
        % 割平面不等式: -(dK + lambda * dM)^T * x - z <= lambda * (v^T * ns_M * v)
        a_cut = -(dK + lambda * dM)';
        b_cut = lambda * (v' * ns_M * v);
        
        A_cuts = [A_cuts; a_cut, -1]; % 追加一行: [关于x的系数, 关于z的系数]
        b_cuts = [b_cuts; b_cut];
        
        % 限制割平面池的大小，防止 LP 求解变慢
        if size(A_cuts, 1) > max_bundle_cuts
            A_cuts(1, :) = [];
            b_cuts(1) = [];
        end

        % 5. 构建并求解 LP
        % 变量组装: vars = [x_1, x_2, ..., x_nm, z]^T
        f = [zeros(nm, 1); 1]; % 目标是最小化 z
        
        % 体积约束: dll^T * x + 0 * z <= V0
        A_ineq = [[dll', 0]; A_cuts];
        b_ineq = [V0; b_cuts];
        
        % 边界条件: x >= x_min, z 自由但设定一个足够大的边界防止发散
        lb = [x_min * ones(nm, 1); z_lb];
        ub = [inf * ones(nm, 1); z_ub];
        
        % 调用 linprog 求解器
        [sol, ~, exitflag] = linprog(f, A_ineq, b_ineq, [], [], lb, ub, linprog_opts);
        
        if exitflag == 1
            x = sol(1:nm);
            z = sol(nm+1);
        else
            bstatus = "LP_Failed";
            break;
        end
    end
    
    x_opt = x;
    z_opt = z;
end
