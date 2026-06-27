clear
% clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Flag.save = 0;
nx = 8;
ny = 8;
maxiter = 40;
% opt_value = -22.2214971322; % nx=6, ny=3, sqrt(2), top right
%opt_value = -51.4026958936; % nx=4, ny=4, Inf, bottom center
% opt_value = -14.5957287796; % nx=8, ny=4, Inf, top right
% opt_value = -3.6909768814; % nx=8, ny=8, Inf, top right
% opt_value = -15.4046034440; % nx=8, ny=6, Inf, bottom center
opt_value = -15.3984958306; % nx=8, ny=8, Inf, bottom center
%opt_value = -27.3905759677; % nx=6, ny=6, Inf, bottom center
lower = opt_value + opt_value*1e-3;
upper = opt_value - opt_value*1e-3;
%% Material constants
V0 = 0.1;
x_min = 1e-8;
scale0 = 1;
eeee = 20000*scale0;
rho = 7.86*10^(-4)*scale0;
dmm0 = 1*scale0;
epsilon = 1e-8;
% [dll,matH,coord_x,ir,irr,ird,Idx] = member(nx,ny,sqrt(2));
[dll,matH,coord_x,ir,irr,ird,Idx] = member(nx,ny,Inf);
% load("member_data44.mat")
%
nk = size(coord_x,1);
nd = size(matH,1); % Number of degree of freedom (size of matrix n)
nm = size(matH,2); % Number of optimization variables m
if Flag.save == 1
    print('-depsc2', '-vector', 'initial_truss');
end
dmm = zeros(nd,1);
% dmm((nd-1):nd,1) = dmm0 * ones(2,1); % top right (original)
% dmm((nd-nx-1):(nd-nx),1) = dmm0 * ones(2,1); % top center
dmm((nx-1):(nx),1) = dmm0 * ones(2,1); % bottom center
ns_M = sparse(diag(dmm));
sqrtK = matH * sparse( diag(sqrt(eeee ./ dll)) );
dm = elem_matrix(coord_x,irr,dll,rho,nm);
strm_M = mk_matrix(dm,ir,nm,nd);
%% Optimization
fprintf(' =========================================================== \n');
fprintf(' ==  Optim. under volume constraint >>> \n');
fprintf('     #members = %g ;  #DOF = %g\n', nm, nd );
iter = 0;
lambdaL = -1000;
lambdaU = 0;
% lambdaL = 0;
% lambdaU = 1;
opt = 1;
t_total = 0;
tic
%%%% Bisection %%%%
% while iter < maxiter
while opt > 1e-4
    tic
    lambda = (lambdaL+lambdaU)/2;
    opt = abs(lambda-opt_value)/abs(opt_value);
    iter = iter+1;
    %%%% SDP subproblem %%%%
    cvx_begin sdp quiet
        cvx_solver sdpt3 % sdpt3 works better than sedumi
        cvx_precision best
%         variable x(nm) nonnegative;
        variables x(nm) z(1); %
        minimize(z); %
        subject to %
        mat_K = sqrtK * sparse(diag(x)) * sqrtK';
        mat_M = reshape(strm_M*x,nd,nd);
%         lambda * mat_K - (mat_M + ns_M) + z * eye(nd) >= 0;
        mat_K + lambda * (mat_M + ns_M) + z * eye(nd) >= 0;
%         z >= 0;
        dll' * x <= V0;
        x >= x_min;
    cvx_end
    %%%%
%     if sum(isnan(x)) >= 1
%         lambdaU = lambda;
%         flag = 1;
%     else
%         lambdaL = lambda;
%         x_new = x;
%         flag = 0;
%         [~] = draw_cs(coord_x,irr,x_new,1);
%     end
%     fprintf('Iter.:%4i Obj.:%15.10f NaN:%1i \n',iter,lambda,flag); 
    %%%%
    if z >= 0
        lambdaL = lambda;
        flag = 0;
    else
        lambdaU = lambda;
        flag = 1;
    end
    t_iter = toc;
    t_total = t_total+t_iter;
%     [~] = draw_cs_freq(coord_x,irr,x,1,nx,ny);
    fprintf('Iter.:%4i Obj.:%15.10f NaN:%1i Time:%8.2f \n',iter,lambda,flag,t_total); 
    %%%%
end
time = toc; %
[~] = draw_cs_freq(coord_x,irr,x,1,nx,ny);
%% Post-process
% fun = @(mult) dll'*max(x_new-mult*dll,x_min)-V0;
% mult_int = [min(x_new-V0)/max(dll),max(x_new)/min(dll)];
% mult = fzero(fun,mult_int);
% x = max(x_new-mult*dll,x_min);
% mat_K = sqrtK*sparse(diag(x))*sqrtK';
% mat_M = reshape(strm_M*x,nd,nd);
% [V,lambda_post,flag] = eigs(-mat_K,mat_M+ns_M,1,"smallestabs");
% [~] = draw_cs(coord_x,irr,x_new,1);
% [V,lambda01,flag] = eigs(-mat_K,mat_M+ns_M,3,"smallestabs");
% mat_K = sqrtK*sparse(diag(x_new))*sqrtK';
% mat_M = reshape(strm_M*x_new,nd,nd);
% [V,lambda02,flag] = eigs(-mat_K,mat_M+ns_M,3,"smallestabs");
%% Output
% fprintf(' lambda = %15.10f [s] \n',-lambda_post); %
fprintf(' time = %6.2f [s] \n',time); %
% [~] = draw_cs(coord_x,irr,x_new,1);
fprintf(' =========================================================== \n');