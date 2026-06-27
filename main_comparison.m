clearvars
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 五大算法全景对比：S-APG, S-PG, Subgrad, CVX(内点法) vs. ManiSDP
% 横轴：问题规模 (4x4, 6x6, 8x8)
% 纵轴：计算时间 (评估规模扩展性) & 最终精度 (评估求解质量)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 预设测试规模与对应的绝对理论最优解 (使用 bottom center, Inf)
scale_list = [4, 6, 8]; 
opt_values = [51.4026958936, 27.3905759677, 15.3984958306]; 

num_scales = length(scale_list);
num_algs = 5; % 1:S-APG, 2:S-PG, 3:Subgrad, 4:CVX(IPM), 5:ManiSDP
alg_names = ["S-APG", "S-PG", "Subgrad", "CVX (IPM)", "ManiSDP"];

final_errors = zeros(num_scales, num_algs);
final_times  = zeros(num_scales, num_algs);

% 传统算法迭代预算
maxiter_trad = 150; 

fprintf('===========================================================\n');
fprintf('开始 [内点法 vs 第一阶 vs ManiSDP] 自动化规模对比...\n');
fprintf('===========================================================\n');

if exist('manopt', 'dir'), addpath(genpath(fullfile(pwd, 'manopt', 'manopt'))); end
warning('off', 'manopt:getHessian:approx');

for s_idx = 1:num_scales
    nx = scale_list(s_idx); ny = nx;
    true_opt = opt_values(s_idx);
    fprintf('\n>>> 正在测试规模: %dx%d (理论基频: %.4f)\n', nx, ny, true_opt);
    
    %% 1. 物理参数与装配
    V0 = 0.1; x_min = 1e-8; scale0 = 1; len_scale = 1;
    eeee = 20000 * scale0; rho = 7.86 * 10^(-4) * scale0; dmm0 = 1 * scale0; epsilon = 1e-8;

    [dll, matH, coord_x, ir, irr, ird, Idx] = member(nx, ny, Inf);
    nd = size(matH, 1); nm = size(matH, 2);

    dmm = zeros(nd, 1); dmm((nx-1):(nx), 1) = dmm0 * ones(2, 1);
    ns_M = sparse(diag(dmm));

    sqrtK = matH * sparse(diag(sqrt(eeee ./ dll)));
    dm = elem_matrix(coord_x, irr, dll, rho, nm);
    strm_M = mk_matrix(dm, ir, nm, nd);
    x0 = V0 * ones(nm, 1) / sum(dll); 
    
    ind_0 = 1:nd^2; ind_1 = ones(nd^2,1); ind_I1 = repmat((1:nd)',nd,1); ind_I2 = zeros(nd^2,1);
    for i = 1:nd, ind_I2((1:nd)+(i-1)*nd) = i*ones(nd,1); end
    I1 = sparse(ind_0,ind_I1,ind_1); I2 = sparse(ind_0,ind_I2,ind_1);

    %% 2. 算法 1-3: 第一阶方法 (S-APG, S-PG, Subgrad)
    alg_to_run = [1, 3, 4]; mu0 = 0.01; L0 = 1000; L0_s = 1000 * 5000; L0_a = 100;
    for idx = 1:length(alg_to_run)
        i = alg_to_run(idx); t_start = tic;
        x = x0; z = x0; a = 0; final_lambda = 0;
        for iter = 1:maxiter_trad
            mat_K = sqrtK * sparse(diag(x)) * sqrtK'; mat_M = reshape(strm_M * x, nd, nd);
            K = full(mat_K); K = K + K'/2; M = full(mat_M + ns_M); M = M + M'/2;
            [~, mat_lambda] = eigs(mat_M + ns_M, mat_K + epsilon * eye(nd), 3);
            final_lambda = 1 / mat_lambda(1,1);
            if i == 1 % S-APG
                a = (1+sqrt(4*a^2+1))/2; y = (1-1/a)*x+(1/a)*z; mu = mu0*(iter)^(-1);
                mat_K_y = sqrtK*sparse(diag(y))*sqrtK'; mat_M_y = reshape(strm_M*y,nd,nd);
                K_y = full(mat_K_y); K_y = K_y+K_y'/2; M_y = full(mat_M_y+ns_M); M_y = M_y+M_y'/2;
                [V_y,lambda_y] = eig(M_y,K_y+epsilon*eye(nd),'vector');
                [lambda_y,ind_y] = sort(lambda_y,'descend'); V_y = V_y(:,ind_y);
                VKV_y = diag(V_y'*(K_y+epsilon*eye(nd))*V_y);
                dK_all_y = (sqrtK'*V_y).^2; V2_y = (I1*V_y).*(I2*V_y); dM_all_y = strm_M'*V2_y;
                df_all_y = dM_all_y-dK_all_y*diag(lambda_y); df_all_y_n = df_all_y./VKV_y';
                exp_dif_y = exp((lambda_y-lambda_y(1))/mu); df_s_y = df_all_y_n*exp_dif_y/sum(exp_dif_y);
                L = L0_a*mu0; z = (z-(a*mu/L)*df_s_y);
                fun = @(mult) dll'*max(z-mult*dll,x_min)-V0; z = max(z-fzero(fun, [min(z-V0)/max(dll),max(z)/min(dll)])*dll,x_min); 
                x_new = (1-1/a)*x+(1/a)*z;
            elseif i == 3 % S-PG
                [V_eig,lambda_eig] = eig(-K,M,'vector'); [lambda_eig,ind_eig] = sort(lambda_eig,'descend'); V_eig = V_eig(:,ind_eig);
                VMV = diag(V_eig'*M*V_eig); mu = mu0*(iter)^(-0.5);
                dK_all = -(sqrtK'*V_eig).^2; V2 = (I1*V_eig).*(I2*V_eig); dM_all = strm_M'*V2;
                df_all = dK_all-dM_all*diag(lambda_eig); df_all_n = df_all./VMV';
                exp_dif = exp((lambda_eig-lambda_eig(1))/mu); df_s = df_all_n*exp_dif/sum(exp_dif);
                L = L0_s*mu0; x_new = x-df_s*mu/L;
                fun = @(mult) dll'*max(x_new-mult*dll,x_min)-V0; x_new = max(x_new-fzero(fun, [min(x_new-V0)/max(dll),max(x_new)/min(dll)])*dll,x_min);
            elseif i == 4 % Subgrad
                [V_eigs,mat_lambda_s] = eigs(mat_M+ns_M,mat_K+epsilon*eye(nd),1);
                lambda_s = mat_lambda_s(1,1); dK = -(sqrtK'*V_eigs(:,1)).^2;
                dM = strm_M'*reshape(V_eigs(:,1)*V_eigs(:,1)',nd*nd,1); df = dK-lambda_s*dM;
                df_n = df/norm(df); L = L0*(iter)^(0.5); x_new = x-df_n/L;
                fun = @(mult) dll'*max(x_new-mult*dll,x_min)-V0; x_new = max(x_new-fzero(fun, [min(x_new-V0)/max(dll),max(x_new)/min(dll)])*dll,x_min);
            end
            x = x_new;
        end
        final_times(s_idx, idx) = toc(t_start);
        final_errors(s_idx, idx) = abs(final_lambda - true_opt) / true_opt;
        fprintf('  - %-10s | 耗时: %6.2f s | 最终误差: %.2e\n', alg_names(idx), final_times(s_idx, idx), final_errors(s_idx, idx));
    end

    %% 3. 算法 4: CVX (内点法) - 二分法外层
    t_start_cvx = tic;
    lambdaL_cvx = -true_opt * 1.2; lambdaU_cvx = 0;
    opt_cvx = 1; iter_cvx = 0; final_lambda_cvx = 0;
    while (lambdaU_cvx - lambdaL_cvx) > 1e-4 && opt_cvx > 1e-4 && iter_cvx < 40
        iter_cvx = iter_cvx + 1;
        lambda = 0.5 * (lambdaL_cvx + lambdaU_cvx);
        opt_cvx = abs(abs(lambda) - true_opt) / true_opt;
        
        cvx_begin sdp quiet
            cvx_solver sdpt3
            variables x_cvx(nm) z_cvx(1);
            minimize(z_cvx);
            subject to
            mat_K = sqrtK * sparse(diag(x_cvx)) * sqrtK';
            mat_M = reshape(strm_M*x_cvx,nd,nd);
            mat_K + lambda * (mat_M + ns_M) + z_cvx * eye(nd) >= 0;
            dll' * x_cvx <= V0;
            x_cvx >= x_min;
        cvx_end
        
        if z_cvx >= 0, lambdaL_cvx = lambda; else, lambdaU_cvx = lambda; end
        final_lambda_cvx = abs(lambdaL_cvx);
    end
    final_times(s_idx, 4) = toc(t_start_cvx);
    final_errors(s_idx, 4) = abs(final_lambda_cvx - true_opt) / true_opt;
    fprintf('  - %-10s | 耗时: %6.2f s | 最终误差: %.2e\n', alg_names(4), final_times(s_idx, 4), final_errors(s_idx, 4));

    %% 4. 算法 5: ManiSDP (流形投影 + 动量)
    t_start_mani = tic;
    lambdaL_mani = -true_opt * 1.2; lambdaU_mani = 0;
    rank_p = min(12, max(4, round(nd / 15))); 
    max_alt = 250; step0_x = 0.8; mu_smooth = 1; beta_mom = 0.90; 
    manifold = stiefelfactory(nd, rank_p); probY.M = manifold; optsY.verbosity = 0; optsY.maxiter = 15; optsY.maxtime = 1.5;

    x_warm = x0; Y_opt = manifold.rand(); iter_mani = 0; final_lambda_mani = 0;
    while (lambdaU_mani - lambdaL_mani) > 1e-4 && iter_mani < 40
        iter_mani = iter_mani + 1;
        lambda = 0.5 * (lambdaL_mani + lambdaU_mani);
        x = x_warm; min_eig = -inf; v_x = zeros(nm, 1); 
        for alt = 1:max_alt
            mat_K = sqrtK * spdiags(x, 0, nm, nm) * sqrtK'; mat_M = reshape(strm_M * x, nd, nd) + ns_M;
            S_mat = mat_K + lambda * mat_M; S_mat = 0.5 * (S_mat + S_mat');
            probY.cost = @(Y) real(trace(Y' * S_mat * Y)); probY.egrad = @(Y) 2 * S_mat * Y;
            [Y_opt, ~, ~, ~] = conjugategradient(probY, Y_opt, optsY);
            H_sub = 0.5 * (Y_opt' * S_mat * Y_opt + (Y_opt' * S_mat * Y_opt)');
            [V_sub, D_sub] = eig(H_sub); eig_sub = real(diag(D_sub)); min_eig = min(eig_sub);
            if min_eig >= -1e-5, break; end
            V_full = Y_opt * V_sub; mu_alt = mu_smooth / sqrt(alt);
            weights = exp(-(eig_sub - min_eig) / (mu_alt + 1e-12)); weights = weights / sum(weights);
            grad_x = zeros(nm, 1);
            for j = 1:rank_p
                vj = V_full(:, j); dK_j = (sqrtK' * vj).^2; dM_j = strm_M' * reshape(vj * vj', nd^2, 1);
                grad_x = grad_x + weights(j) * (dK_j + lambda * dM_j);
            end
            v_x = beta_mom * v_x + (1 - beta_mom) * (grad_x / (norm(grad_x) + 1e-12));
            step = step0_x / (1 + 0.05 * alt); x_new = max(real(x + step * v_x * V0), x_min);
            fun = @(mult) dll'*max(x_new-mult*dll,x_min)-V0; x = max(x_new-fzero(fun, [min(x_new-V0)/max(dll)-1, max(x_new)/min(dll)+1])*dll,x_min);
        end
        x_warm = x;
        if -min_eig >= 0, lambdaL_mani = lambda; else, lambdaU_mani = lambda; end
        final_lambda_mani = abs(lambdaL_mani);
    end
    final_times(s_idx, 5) = toc(t_start_mani);
    final_errors(s_idx, 5) = abs(final_lambda_mani - true_opt) / true_opt;
    fprintf('  - %-10s | 耗时: %6.2f s | 最终误差: %.2e\n', alg_names(5), final_times(s_idx, 5), final_errors(s_idx, 5));
end

%% 4. 绘制双子图：规模 vs 耗时 & 规模 vs 精度
fprintf('\n===========================================================\n');
final_errors = max(final_errors, 1e-12); % 防止 log(0)
color_map = ["#0072BD", "#EDB120", "#7E2F8E", "#77AC30", "#D95319"]; 
line_styles = {'--o', '--s', '--^', '-d', '-p'};

figure('Name', 'Comprehensive Comparison', 'Position', [100, 100, 1400, 550]);

% 子图 1：耗时对比
subplot(1, 2, 1);
for idx = 1:num_algs
    line_w = 2; marker_s = 8;
    if idx == 4 || idx == 5, line_w = 3; marker_s = 10; end % 加粗 CVX 和 ManiSDP
    semilogy(scale_list, final_times(:, idx), line_styles{idx}, 'LineWidth', line_w, ...
        'MarkerSize', marker_s, 'Color', color_map(idx), 'MarkerFaceColor', color_map(idx)); hold on;
end
grid on; box on; set(gca, 'XTick', scale_list, 'XTickLabel', arrayfun(@(x) sprintf('%dx%d', x, x), scale_list, 'UniformOutput', false));
xlabel('Grid Size', 'FontSize', 14); ylabel('Wall-clock Time (Seconds) [Log Scale]', 'FontSize', 14);
title('Computational Scalability (Time vs Scale)', 'FontSize', 16);
legend(alg_names, 'Location', 'northwest', 'FontSize', 12);

% 子图 2：精度对比
subplot(1, 2, 2);
for idx = 1:num_algs
    line_w = 2; marker_s = 8;
    if idx == 4 || idx == 5, line_w = 3; marker_s = 10; end
    semilogy(scale_list, final_errors(:, idx), line_styles{idx}, 'LineWidth', line_w, ...
        'MarkerSize', marker_s, 'Color', color_map(idx), 'MarkerFaceColor', color_map(idx)); hold on;
end
grid on; box on; set(gca, 'XTick', scale_list, 'XTickLabel', arrayfun(@(x) sprintf('%dx%d', x, x), scale_list, 'UniformOutput', false));
xlabel('Grid Size', 'FontSize', 14); ylabel('Final Relative Error [Log Scale]', 'FontSize', 14);
title('Solution Precision (Error vs Scale)', 'FontSize', 16);
ylim([1e-6, 1e-0]);

sgtitle('Comprehensive Performance: First-order vs. IPM vs. ManiSDP', 'FontSize', 18, 'FontWeight', 'bold');