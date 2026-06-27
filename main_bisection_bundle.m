function [x_center, z_best, best_min_eig, bundle_size, status] = solve_sdp_bundle( ...
    lambda, x0, dll, V0, x_min, sqrtK, strm_M, ns_M, nd, nm, ...
    tol_psd, tol_z, max_iter, max_bundle_size, z_lower, z_upper, linprog_opts)
    
    x = x0;
    x_center = x0; % 信赖域中心
    
    % 先评估初始点
    S_center = assemble_s_matrix(x_center, 0, lambda, sqrtK, strm_M, ns_M, nd);
    [best_min_eig, ~] = smallest_eigpair(S_center);
    
    % 快速通道：如果初始点已经满足，直接跳出
    if best_min_eig >= -tol_psd
        z_best = -best_min_eig;
        bundle_size = 0;
        status = "psd-feasible";
        return;
    end
    
    delta = 0.05 * V0; % 初始信赖域半径
    cut_C = []; cut_alpha = []; cut_const = [];
    status = "maxcut";
    z_lp = 0;
    
    for k = 1:max_iter
        % 评估当前点的真实矩阵 (不带 z)
        S_true = assemble_s_matrix(x, 0, lambda, sqrtK, strm_M, ns_M, nd);
        [true_min_eig, v] = smallest_eigpair(S_true);
        
        if true_min_eig >= -tol_psd
            status = "psd-feasible";
            best_min_eig = true_min_eig;
            x_center = x;
            z_lp = -true_min_eig;
            break;
        end
        
        % 添加当前点的割线
        [c, alpha, beta] = make_cut(v, lambda, sqrtK, strm_M, ns_M, nd);
        
        % 【核心稳定器 1】割线归一化，拯救内点法的病态矩阵
        cut_norm = norm([c; alpha; beta]) + 1e-12;
        c = c / cut_norm; alpha = alpha / cut_norm; beta = beta / cut_norm;
        
        cut_C = [cut_C; c']; cut_alpha = [cut_alpha; alpha]; cut_const = [cut_const; beta];
        if size(cut_C, 1) > max_bundle_size
            cut_C = cut_C(end-max_bundle_size+1:end, :);
            cut_alpha = cut_alpha(end-max_bundle_size+1:end);
            cut_const = cut_const(end-max_bundle_size+1:end);
        end
        
        % 【核心稳定器 2】信赖域动态更新策略 (Serious Step vs Null Step)
        if true_min_eig > best_min_eig + 1e-5
            best_min_eig = true_min_eig; % 确有进步，更新中心点
            x_center = x;
            delta = min(0.2 * V0, delta * 1.2); % 扩大搜索范围
        else
            delta = max(1e-4 * V0, delta * 0.7); % 原地踏步，缩小搜索范围逼近真实解
        end
        
        % 设定信赖域边界 (限制 LP 的横跳范围)
        lb_x = max(x_min, x_center - delta);
        ub_x = x_center + delta; 
        
        % 求解带信赖域的线性规划
        [x_lp, z_lp, exitflag] = solve_cut_lp(cut_C, cut_alpha, cut_const, ...
            dll, V0, lb_x, ub_x, z_lower, z_upper, linprog_opts);
            
        if exitflag <= 0
            status = "lp-failed";
            break;
        end
        
        x = x_lp;
        
        if z_lp > tol_z
            status = "infeasible-z";
            break;
        end
    end
    bundle_size = size(cut_C, 1);
    z_best = z_lp;
end

function [x, z, exitflag] = solve_cut_lp(cut_C, cut_alpha, cut_const, ...
    dll, V0, lb_x, ub_x, z_lower, z_upper, opts)
    
    f = [zeros(length(lb_x), 1); 1];
    A = [dll', 0; -cut_C, -cut_alpha];
    b = [V0; cut_const];
    lb = [lb_x; z_lower];
    ub = [ub_x; z_upper];
    
    [sol, ~, exitflag] = linprog(f, A, b, [], [], lb, ub, opts);
    
    if exitflag <= 0 || isempty(sol) || any(isnan(sol)) || any(isinf(sol))
        x = [];
        z = inf;
        exitflag = -1;
    else
        x = sol(1:end-1);
        z = sol(end);
    end
end

function [min_eig, v] = smallest_eigpair(S)
    if any(isnan(S(:))) || any(isinf(S(:)))
        min_eig = -inf;
        v = ones(size(S, 1), 1) / sqrt(size(S, 1));
        return;
    end
    S = full(0.5 * (S + S'));
    [V, D] = eig(S, 'vector');
    [min_eig, idx] = min(real(D));
    v = V(:, idx);
    v = real(v);
    if norm(v) > 1e-12
        v = v / norm(v);
    else
        v = ones(size(S, 1), 1) / sqrt(size(S, 1));
    end
end

function [c, alpha, beta] = make_cut(v, lambda, sqrtK, strm_M, ns_M, nd)
    cK = (sqrtK' * v).^2;
    vv = reshape(v * v', nd^2, 1);
    cM = strm_M' * vv;
    c = cK + lambda * cM;
    alpha = v' * v;
    beta = lambda * (v' * ns_M * v);
end

function S = assemble_s_matrix(x, z, lambda, sqrtK, strm_M, ns_M, nd)
    mat_K = sqrtK * sparse(diag(x)) * sqrtK';
    mat_M = reshape(strm_M * x, nd, nd);
    S = mat_K + lambda * (mat_M + ns_M) + z * speye(nd);
    S = 0.5 * (S + S');
end

function x = project_volume(x_raw, dll, V0, x_min)
    x_raw = max(real(x_raw), x_min);
    lo = min((x_raw - V0) ./ dll) - 1;
    hi = max(x_raw ./ dll) + 1;
    for k = 1:80
        mid = 0.5 * (lo + hi);
        x_mid = max(x_raw - mid * dll, x_min);
        if dll' * x_mid > V0
            lo = mid;
        else
            hi = mid;
        end
    end
    x = max(x_raw - hi * dll, x_min);
end