

clearvars;
close all;
clc;

%% ManiSDP / Burer-Monteiro ALM for eigenfrequency maximization
% This script keeps the structural model used in main_bisection.m and
% replaces the CVX SDP feasibility test by a low-rank PSD factorization.

rng(1);
warning('off', 'manopt:getHessian:approx');

if exist('manopt', 'dir')
    addpath(genpath(fullfile(pwd, 'manopt', 'manopt')));
end

%% Structural parameters
Flag.save = 0;
nx = 4;
ny = 4;

V0 = 0.1;
x_min = 1e-8;
scale0 = 1;
eeee = 20000 * scale0;
rho = 7.86e-4 * scale0;
dmm0 = 1 * scale0;
epsilon = 1e-8;

[dll, matH, coord_x, ir, irr, ird, Idx] = member(nx, ny, Inf);
nd = size(matH, 1);
nm = size(matH, 2);

if Flag.save == 1
    print('-depsc2', '-vector', 'initial_truss');
end

dmm = zeros(nd, 1);
dmm((nx-1):nx, 1) = dmm0 * ones(2, 1);  % bottom center non-structural mass
ns_M = sparse(diag(dmm));

sqrtK = matH * sparse(diag(sqrt(eeee ./ dll)));
dm = elem_matrix(coord_x, irr, dll, rho, nm);
strm_M = mk_matrix(dm, ir, nm, nd);

%% Algorithmic parameters
lambdaL = -1000;
lambdaU = 0;
tol_bisect = 1e-3;
maxiter_bisect = 25;

rank_r = min(20, max(8, ceil((sqrt(8 * nd + 1) - 1) / 2) + 2));
maxiter_alm = 120;
maxiter_x = 20;
maxiter_linesearch = 12;
tol_psd_res = 5e-3;
tol_psd_eig = 1e-7;
rho0 = 0.5;
rho_max = 1e4;
step0_x = 0.2;
step_shrink = 0.5;
armijo_c = 1e-4;

fprintf(' =========================================================== \n');
fprintf(' ==  ManiSDP / BM-ALM under equality volume constraint >>> \n');
fprintf('     #members = %g ; #DOF = %g ; rank = %g\n', nm, nd, rank_r);

if exist('symfixedrankYYfactory', 'file') ~= 2 || exist('trustregions', 'file') ~= 2
    error(['Manopt is required for main_manisdp.m. ', ...
        'Please make sure the manopt folder is available and on the MATLAB path.']);
end

manifold = symfixedrankYYfactory(nd, rank_r);
problem.M = manifold;
options.verbosity = 0;
options.maxiter = 8;
options.maxtime = 0.8;
options.minstepsize = 1e-12;

%% Initial design
x0 = ones(nm, 1) + 0.05 * randn(nm, 1);
x0 = project_volume(x0, dll, V0, x_min);
x = x0;
x_best = x;
V_opt = manifold.rand();

t_total = 0;
iter = 0;
history = struct('lambda', [], 'residual', [], 'feasible', [], 'time', []);

while abs(lambdaU - lambdaL) > tol_bisect && iter < maxiter_bisect
    t_iter_start = tic;
    iter = iter + 1;
    lambda = 0.5 * (lambdaL + lambdaU);

    U = zeros(nd, nd);
    rho_pen = rho0;
    old_residual = Inf;
    x_trial = x_best;
    z_trial = 0;

    for alm_iter = 1:maxiter_alm
        [mat_K, mat_M, S_base] = assemble_sdp_matrix(x_trial, lambda, sqrtK, strm_M, ns_M, nd);
        S_mat = S_base + z_trial * speye(nd);

        A_tgt = S_mat + U / rho_pen;
        A_tgt = 0.5 * (A_tgt + A_tgt');

        problem.cost = @(V) 0.5 * norm(V * V' - A_tgt, 'fro')^2;
        problem.egrad = @(V) 2 * (V * V' - A_tgt) * V;
        problem.ehess = @(V, dV) 2 * ((dV * V' + V * dV') * V + (V * V' - A_tgt) * dV);
        [V_opt, ~, ~, ~] = trustregions(problem, V_opt, options);

        A_tgt_x = V_opt * V_opt' - U / rho_pen;
        for inner = 1:maxiter_x
            [~, ~, S_base] = assemble_sdp_matrix(x_trial, lambda, sqrtK, strm_M, ns_M, nd);
            S_mat = S_base + z_trial * speye(nd);
            W = S_mat - A_tgt_x;
            cost_xz = z_trial / rho_pen + 0.5 * norm(W, 'fro')^2;

            grad_K = sum((W * sqrtK) .* sqrtK, 1)';
            grad_M = strm_M' * W(:);
            grad_x = grad_K + lambda * grad_M;
            grad_z = trace(W) + 1 / rho_pen;

            grad_norm = norm([grad_x; grad_z]);
            if grad_norm < 1e-12
                break;
            end

            step_size = min(step0_x, (0.08 * V0) / grad_norm);
            accepted = false;
            for ls_iter = 1:maxiter_linesearch
                x_candidate = project_volume(x_trial - step_size * grad_x, dll, V0, x_min);
                z_candidate = z_trial - step_size * grad_z;
                [~, ~, S_candidate] = assemble_sdp_matrix(x_candidate, lambda, sqrtK, strm_M, ns_M, nd);
                S_candidate = S_candidate + z_candidate * speye(nd);
                W_candidate = S_candidate - A_tgt_x;
                cost_candidate = z_candidate / rho_pen + 0.5 * norm(W_candidate, 'fro')^2;

                if cost_candidate <= cost_xz - armijo_c * step_size * grad_norm^2 || cost_candidate < cost_xz
                    x_trial = x_candidate;
                    z_trial = z_candidate;
                    accepted = true;
                    break;
                end

                step_size = step_shrink * step_size;
            end

            if ~accepted
                break;
            end
        end

        [~, ~, S_base] = assemble_sdp_matrix(x_trial, lambda, sqrtK, strm_M, ns_M, nd);
        S_mat = S_base + z_trial * speye(nd);
        residual_mat = S_mat - V_opt * V_opt';
        residual = norm(residual_mat, 'fro') / max(1, norm(S_mat, 'fro'));
        min_sdp_eig = min(real(eig(full(S_mat))));

        if residual < tol_psd_res && min_sdp_eig >= -tol_psd_eig
            break;
        end

        U = U + rho_pen * residual_mat;
        if residual > 0.7 * old_residual
            rho_pen = min(1.5 * rho_pen, rho_max);
        end
        old_residual = residual;
    end

    feasible = residual < tol_psd_res && min_sdp_eig >= -tol_psd_eig && z_trial < 0;
    if feasible
        lambdaU = lambda;
        x_best = x_trial;
        x = x_trial;
    else
        lambdaL = lambda;
        x = 0.8 * x_best + 0.2 * x0;
        x = project_volume(x, dll, V0, x_min);
    end

    t_total = t_total + toc(t_iter_start);
    history.lambda(iter, 1) = lambda;
    history.residual(iter, 1) = residual;
    history.feasible(iter, 1) = feasible;
    history.time(iter, 1) = t_total;
    history.min_sdp_eig(iter, 1) = min_sdp_eig;
    history.z(iter, 1) = z_trial;

    fprintf('Iter.:%4i Obj.:%15.10f Feas:%1i z:%10.3e ALM:%4i Res:%9.2e MinEig:%9.2e Time:%8.2f \n', ...
        iter, lambda, feasible, z_trial, alm_iter, residual, min_sdp_eig, t_total);
end

%% Post-process
[lambda_freq, lambda_all] = fundamental_frequency_value(x_best, sqrtK, strm_M, ns_M, nd, epsilon);

fprintf(' =========================================================== \n');
fprintf(' ManiSDP finished. \n');
fprintf(' SDP bisection value     = %15.10f \n', -lambdaU);
fprintf(' eigenvalue check Eig1   = %15.10f \n', lambda_freq);
fprintf(' volume                  = %15.10f / %15.10f \n', dll' * x_best, V0);
fprintf(' total time              = %8.2f [s] \n', t_total);
fprintf(' =========================================================== \n');

figure(1);
[~] = draw_cs_freq(coord_x, irr, x_best, 1, nx, ny);
title(sprintf('ManiSDP Final Topology (%dx%d), Eig1 = %.4f', nx, ny, lambda_freq));

figure(2);
semilogy(history.residual, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Bisection iteration');
ylabel('relative PSD residual');
title('ManiSDP feasibility residual');

save('result_manisdp.mat', 'x_best', 'lambdaU', 'lambda_freq', 'lambda_all', ...
    'history', 'nx', 'ny', 'V0', 'x_min', 'rank_r');

%% Local functions
function x = project_volume(x_raw, dll, V0, x_min)
    x_raw = max(real(x_raw), x_min);

    min_volume = x_min * sum(dll);
    if min_volume > V0
        error('The lower bound x_min is infeasible: x_min * sum(dll) > V0.');
    end

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

function [mat_K, mat_M, S_mat] = assemble_sdp_matrix(x, lambda, sqrtK, strm_M, ns_M, nd)
    mat_K = sqrtK * sparse(diag(x)) * sqrtK';
    mat_M = reshape(strm_M * x, nd, nd);
    S_mat = mat_K + lambda * (mat_M + ns_M);
    S_mat = 0.5 * (S_mat + S_mat');
end

function [lambda_1, lambda_all] = fundamental_frequency_value(x, sqrtK, strm_M, ns_M, nd, epsilon)
    mat_K = sqrtK * sparse(diag(x)) * sqrtK';
    mat_M = reshape(strm_M * x, nd, nd);
    K = full(0.5 * (mat_K + mat_K') + epsilon * eye(nd));
    M = full(0.5 * (mat_M + mat_M') + ns_M);

    try
        lambda_inv = eigs(M, K, 3);
        lambda_inv = real(lambda_inv(:));
        lambda_all = sort(1 ./ lambda_inv, 'ascend');
    catch
        lambda_all = sort(real(eig(K, M)), 'ascend');
    end

    lambda_all = lambda_all(isfinite(lambda_all) & lambda_all > 0);
    lambda_1 = lambda_all(1);
end
