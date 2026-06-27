clearvars
close all;
tic
%% Algorithmic parameters
alg_list = ["S-APG" "Inexact S-APG" "S-PG" "Subgrad"];
% alg_list = ["S-APG" "Inexact S-APG $l=1$" "Inexact S-APG $l=2$" "Inexact S-APG $l=3$"];
% alg_list = ["S-APG" "S-PG ($0.99^k$)" "S-PG ($k^{-1/2}$)" "Subgrad"];
alg = [1 3 4];
n_alg = numel(alg);
alg_name = alg_list(alg);
maxiter = 100;
V0 = 0.1; % 0.011;
x_min = 1e-8;
L0 = 1000;
L0_s = 1000*5000; % S-APG L0=100 (3000iter), S-PG L0=1000 for non-oscillatory
L0_a = 100;
% L0_a = 0.1; % without NS
mu0 = 0.01;
%% Structural parameters
flag = 0;
Flag.save = 0;
nx = 8;
ny = 8;
%%%% Material constants %%%%
scale0 = 1; % aligned with main_bisection.m
len_scale = 1; % m -> 10^(-2)m
eeee = 20000*(len_scale^2)*scale0;
rho = 7.86*10^(-4)*(len_scale^3)*scale0;
dmm0 = 1*scale0;
epsilon = 1e-8;
% dmm0 = 0;
% [dll,matH,coord_x,ir,irr,ird,Idx] = member(nx,ny,sqrt(2));
[dll,matH,coord_x,ir,irr,ird,Idx] = member(nx,ny,Inf);
nk = size(coord_x,1);
nd = size(matH,1); % Number of degree of freedom (size of matrix n)
nm = size(matH,2); % Number of optimization variables m
if Flag.save == 1
    print('-depsc2', '-vector', 'initial_truss');
end
%%%% Non-structural mass %%%%
dmm = zeros(nd,1);
% dmm((nd-1):nd,1) = dmm0 * ones(2,1); % top right (original)
% dmm((nd-nx-1):(nd-nx),1) = dmm0 * ones(2,1); % top center
dmm((nx-1):(nx),1) = dmm0 * ones(2,1); % bottom center
% dmm(1:2*(nx-1),1) = dmm0 * ones(2*(nx-1),1); % bottom all
% dmm((nx-5):(nx+6),1) = dmm0 * ones(12,1); % bottom
ns_M = sparse(diag(dmm));
%%%% Stiffness & mass matrix %%%%
sqrtK = matH * sparse( diag(sqrt(eeee ./ dll)) );
dm = elem_matrix(coord_x,irr,dll,rho,nm);
strm_M = mk_matrix(dm,ir,nm,nd);
% [~] = draw_cs_ini_freq(coord_x, irr, ones(nm,1), nx, ny); % draw initial design
%% Setting for eigenfrequency maximization
ind_0 = 1:nd^2;
ind_1 = ones(nd^2,1);
ind_I1 = repmat((1:nd)',nd,1);
ind_I2 = zeros(nd^2,1);
for i = 1:nd
    ind_I2((1:nd)+(i-1)*nd) = i*ones(nd,1);
end
I1 = sparse(ind_0,ind_I1,ind_1);
I2 = sparse(ind_0,ind_I2,ind_1);
data_lambda = zeros(maxiter+1,numel(alg_list),3);
data_obj = zeros(maxiter,numel(alg_list));
data_ch = zeros(maxiter,numel(alg_list));
t_iter = zeros(maxiter,numel(alg_list));
x0 = V0*ones(nm,1)/sum(dll);
%% Algorithms
fprintf(' =========================================================== \n');
fprintf(' ==  Optim. under volume constraint >>> \n');
fprintf('     #members = %g ;  #DOF = %g\n', nm, nd );
for i = alg
t_alg_start = tic;  % <--- 新增：为当前算法开启独立计时
iter = 0;
x = x0;
%tic
%for i = alg
%iter = 0;
%x = x0;
z = x0;
a = 0;
mu = mu0;
while iter <= maxiter
    iter = iter+1;
    mat_K = sqrtK*sparse(diag(x))*sqrtK';
    mat_M = reshape(strm_M*x,nd,nd);
    K = full(mat_K);
    K = K+K'/2;
    M = full(mat_M+ns_M);
    M = M+M'/2;
%     [V,mat_lambda,flag] = eigs(-mat_K-epsilon*eye(nd),mat_M+ns_M,3,"smallestabs"); % sparse
    [V,mat_lambda,flag] = eigs(mat_M+ns_M,mat_K+epsilon*eye(nd),3); % inverse
    lambda = diag(mat_lambda); % sparse
    lambda = 1./lambda;
%     [V,lambda] = eig(-K,M,'vector');
%     [lambda,ind] = sort(lambda,'descend');
%     V = V(:,ind);
    data_obj(iter,i) = lambda(1);  
    data_lambda(iter,i,:) = lambda(1:3)';
%     data_obj(iter,i) = -lambda(1);  
%     data_lambda(iter,i,:) = -lambda(1:3)';
    tic
    if i == 1
    %% Smoothing accelerated gradient [d'Aspremont2022,Alg20]
    a = (1+sqrt(4*a^2+1))/2;
    y = (1-1/a)*x+(1/a)*z;
    mu = mu0*(iter)^(-1);
    mat_K_y = sqrtK*sparse(diag(y))*sqrtK';
    mat_M_y = reshape(strm_M*y,nd,nd);
    K_y = full(mat_K_y);
    K_y = K_y+K_y'/2;
    M_y = full(mat_M_y+ns_M);
    M_y = M_y+M_y'/2;
%     [V_y,lambda_y] = eig(-K_y-epsilon*eye(nd),M_y,'vector');
    [V_y,lambda_y] = eig(M_y,K_y+epsilon*eye(nd),'vector');
    [lambda_y,ind_y] = sort(lambda_y,'descend');
    V_y = V_y(:,ind_y);
%     VMV_y = diag(V_y'*M_y*V_y);
    VKV_y = diag(V_y'*(K_y+epsilon*eye(nd))*V_y);
    dK_all_y = (sqrtK'*V_y).^2;
    V2_y = (I1*V_y).*(I2*V_y);
    dM_all_y = strm_M'*V2_y;
%     df_all_y = dK_all_y-dM_all_y*diag(lambda_y);
    df_all_y = dM_all_y-dK_all_y*diag(lambda_y);
%     df_all_y_n = df_all_y; % No normalization
%     df_all_y_n = normc(df_all_y); % Normalization
%     df_all_y_n = df_all_y./VMV_y';% Normalization w.r.t. M
    df_all_y_n = df_all_y./VKV_y';% Normalization w.r.t. K
    exp_dif_y = exp((lambda_y-lambda_y(1))/mu);
    df_s_y = df_all_y_n*exp_dif_y/sum(exp_dif_y);
    normdf = norm(df_s_y);
%     df_s_y_n = df_s_y/norm(df_s_y); % Normalization of s-gradient
    df_s_y_n = df_s_y; % No normalization
    L = L0_a*mu0;
    %%%% Adaptive stepsize
%     if iter > 1
% %         L = mu0*norm(df_s_y_n-df_old)/norm(z-z_old);
%         L = (L+norm(df_s_y_n-df_old)^2/abs((df_s_y_n-df_old)'*(z-z_old)))/2;
%     end
%     z_old = z;
%     df_old = df_s_y_n;
    %%%%
    z = (z-(a*mu/L)*df_s_y_n);
    fun = @(mult) dll'*max(z-mult*dll,x_min)-V0;
    mult_int = [min(z-V0)/max(dll),max(z)/min(dll)];
    mult = fzero(fun,mult_int);
    z = max(z-mult*dll,x_min);
    x_new = (1-1/a)*x+(1/a)*z;
    %%%%
    elseif i == 2
    %% Inexact smoothing accelerated gradient
    a = (1+sqrt(4*a^2+1))/2;
    y = (1-1/a)*x+(1/a)*z;
    mu = mu0*(iter)^(-1);
    mat_K_y = sqrtK*sparse(diag(y))*sqrtK';
    mat_M_y = reshape(strm_M*y,nd,nd);
    [V_y,mat_lambda_y] = eigs(-mat_K_y,mat_M_y+ns_M,1,"smallestabs");
    lambda_y = diag(mat_lambda_y);
    [lambda_y,ind_y] = sort(lambda_y,'descend');
    V_y = V_y(:,ind_y);
    VMV_y = diag(V_y'*M_y*V_y);
    dK_all_y = -(sqrtK'*V_y).^2;
    V2_y = (I1*V_y).*(I2*V_y);
    dM_all_y = strm_M'*V2_y;
    df_all_y = dK_all_y-dM_all_y*diag(lambda_y);
%     df_all_y_n = df_all_y; % No normalization
%     df_all_y_n = normc(df_all_y); % Normalization
    df_all_y_n = df_all_y./VMV_y';% Normalization w.r.t. M
    exp_dif_y = exp((lambda_y-lambda_y(1))/mu);
    df_s_y = df_all_y_n*exp_dif_y/sum(exp_dif_y);
    normdf = norm(df_s_y);
%     df_s_y_n = df_s_y/norm(df_s_y); % Normalization of s-gradient
    df_s_y_n = df_s_y; % No normalization
    L = L0_a*mu0;
    z = (z-(a*mu/L)*df_s_y_n);
    fun = @(mult) dll'*max(z-mult*dll,x_min)-V0;
    mult_int = [min(z-V0)/max(dll),max(z)/min(dll)];
    mult = fzero(fun,mult_int);
    z = max(z-mult*dll,x_min);
    x_new = (1-1/a)*x+(1/a)*z;
    %%%% Smoothing method (exponential decay)
%     mat_K = sqrtK*sparse(diag(x))*sqrtK';
%     mat_M = reshape(strm_M*x,nd,nd);
%     K = full(mat_K);
%     K = K+K'/2;
%     M = full(mat_M+ns_M);
%     M = M+M'/2;
%     [V,lambda] = eig(-K,M,'vector'); % dense
%     [lambda,ind] = sort(lambda,'descend');
%     V = V(:,ind);
%     VMV = diag(V'*M*V);
%     mu = mu0*(0.99)^(iter);
%     dK_all = -(sqrtK'*V).^2;
%     V2 = (I1*V).*(I2*V);
%     dM_all = strm_M'*V2;
%     df_all = dK_all-dM_all*diag(lambda);
% %     df_all_n = df_all; % No normalization
% %     df_all_n = normc(df_all); % Normalization
%     df_all_n = df_all./VMV';% Normalization w.r.t. M
%     exp_dif = exp((lambda-lambda(1))/mu);
%     df_s = df_all_n*exp_dif/sum(exp_dif);
%     normdf = norm(df_s);
% %     df_s_n = df_s/norm(df_s); % Normalization of s-gradient
%     df_s_n = df_s; % No normalization
%     L = L0_s*mu0;
%     x_new = x-df_s_n*mu/L;
%     fun = @(mult) dll'*max(x_new-mult*dll,x_min)-V0;
%     mult_int = [min(x_new-V0)/max(dll),max(x_new)/min(dll)];
%     mult = fzero(fun,mult_int);
%     x_new = max(x_new-mult*dll,x_min);
    %%%%
    elseif i == 3
        %% Inexact smoothing accelerated gradient
%     a = (1+sqrt(4*a^2+1))/2;
%     y = (1-1/a)*x+(1/a)*z;
%     mu = mu0*(iter)^(-1);
%     mat_K_y = sqrtK*sparse(diag(y))*sqrtK';
%     mat_M_y = reshape(strm_M*y,nd,nd);
%     [V_y,mat_lambda_y] = eigs(-mat_K_y,mat_M_y+ns_M,2,"smallestabs");
%     lambda_y = diag(mat_lambda_y);
%     [lambda_y,ind_y] = sort(lambda_y,'descend');
%     V_y = V_y(:,ind_y);
%     VMV_y = diag(V_y'*M_y*V_y);
%     dK_all_y = -(sqrtK'*V_y).^2;
%     V2_y = (I1*V_y).*(I2*V_y);
%     dM_all_y = strm_M'*V2_y;
%     df_all_y = dK_all_y-dM_all_y*diag(lambda_y);
% %     df_all_y_n = df_all_y; % No normalization
% %     df_all_y_n = normc(df_all_y); % Normalization
%     df_all_y_n = df_all_y./VMV_y';% Normalization w.r.t. M
%     exp_dif_y = exp((lambda_y-lambda_y(1))/mu);
%     df_s_y = df_all_y_n*exp_dif_y/sum(exp_dif_y);
% %     df_s_y_n = df_s_y/norm(df_s_y); % Normalization of s-gradient
%     df_s_y_n = df_s_y; % No normalization
%     L = L0_a*mu0;
%     z = (z-(a*mu/L)*df_s_y_n);
%     fun = @(mult) dll'*max(z-mult*dll,x_min)-V0;
%     mult_int = [min(z-V0)/max(dll),max(z)/min(dll)];
%     mult = fzero(fun,mult_int);
%     z = max(z-mult*dll,x_min);
%     x_new = (1-1/a)*x+(1/a)*z;
    %%%%  
    %% Smoothing gradient
    mat_K = sqrtK*sparse(diag(x))*sqrtK';
    mat_M = reshape(strm_M*x,nd,nd);
    K = full(mat_K);
    K = K+K'/2;
    M = full(mat_M+ns_M);
    M = M+M'/2;
    [V,lambda] = eig(-K,M,'vector'); % dense
    [lambda,ind] = sort(lambda,'descend');
    V = V(:,ind);
    VMV = diag(V'*M*V);
    mu = mu0*(iter)^(-0.5);
    dK_all = -(sqrtK'*V).^2;
    V2 = (I1*V).*(I2*V);
    dM_all = strm_M'*V2;
    df_all = dK_all-dM_all*diag(lambda);
%     df_all_n = df_all; % No normalization
%     df_all_n = normc(df_all); % Normalization
    df_all_n = df_all./VMV';% Normalization w.r.t. M
    exp_dif = exp((lambda-lambda(1))/mu);
    df_s = df_all_n*exp_dif/sum(exp_dif);
    normdf = norm(df_s);
%     df_s_n = df_s/norm(df_s); % Normalization of s-gradient
    df_s_n = df_s; % No normalization
    L = L0_s*mu0;
    x_new = x-df_s_n*mu/L;
    fun = @(mult) dll'*max(x_new-mult*dll,x_min)-V0;
    mult_int = [min(x_new-V0)/max(dll),max(x_new)/min(dll)];
    mult = fzero(fun,mult_int);
    x_new = max(x_new-mult*dll,x_min);
    %%%%
    elseif i == 4
    %% Inexact smoothing accelerated gradient
%     a = (1+sqrt(4*a^2+1))/2;
%     y = (1-1/a)*x+(1/a)*z;
%     mu = mu0*(iter)^(-1);
%     mat_K_y = sqrtK*sparse(diag(y))*sqrtK';
%     mat_M_y = reshape(strm_M*y,nd,nd);
%     [V_y,mat_lambda_y] = eigs(-mat_K_y,mat_M_y+ns_M,3,"smallestabs");
%     lambda_y = diag(mat_lambda_y);
%     [lambda_y,ind_y] = sort(lambda_y,'descend');
%     V_y = V_y(:,ind_y);
%     VMV_y = diag(V_y'*M_y*V_y);
%     dK_all_y = -(sqrtK'*V_y).^2;
%     V2_y = (I1*V_y).*(I2*V_y);
%     dM_all_y = strm_M'*V2_y;
%     df_all_y = dK_all_y-dM_all_y*diag(lambda_y);
% %     df_all_y_n = df_all_y; % No normalization
% %     df_all_y_n = normc(df_all_y); % Normalization
%     df_all_y_n = df_all_y./VMV_y';% Normalization w.r.t. M
%     exp_dif_y = exp((lambda_y-lambda_y(1))/mu);
%     df_s_y = df_all_y_n*exp_dif_y/sum(exp_dif_y);
% %     df_s_y_n = df_s_y/norm(df_s_y); % Normalization of s-gradient
%     df_s_y_n = df_s_y; % No normalization
%     L = L0_a*mu0;
%     z = (z-(a*mu/L)*df_s_y_n);
%     fun = @(mult) dll'*max(z-mult*dll,x_min)-V0;
%     mult_int = [min(z-V0)/max(dll),max(z)/min(dll)];
%     mult = fzero(fun,mult_int);
%     z = max(z-mult*dll,x_min);
%     x_new = (1-1/a)*x+(1/a)*z;
    %%%%
    %% Subgradient
    mat_K = sqrtK*sparse(diag(x))*sqrtK';
    mat_M = reshape(strm_M*x,nd,nd);
    K = full(mat_K);
    K = K+K'/2;
    M = full(mat_M+ns_M);
    M = M+M'/2;
    [V,mat_lambda_s] = eigs(mat_M+ns_M,mat_K+epsilon*eye(nd),1);
    lambda_s = diag(mat_lambda_s);
    [lambda_s,ind] = sort(lambda_s,'descend');
    V = V(:,ind);
%     vMv1 = V(:,1)'*M*V(:,1);
    dK = -(sqrtK'*V(:,1)).^2;
    dM = strm_M'*reshape(V(:,1)*V(:,1)',nd*nd,1);
    df = dK-lambda_s(1)*dM;
    df_n = df/norm(df); % normalized grad
%     df_n = df/vMv1; % normalized eigenvector w.r.t. M (too large, why?)
    L = L0*(iter)^(0.5);
    x_new = x-df_n/L;
    fun = @(mult) dll'*max(x_new-mult*dll,x_min)-V0;
    mult_int = [min(x_new-V0)/max(dll),max(x_new)/min(dll)];
    mult = fzero(fun,mult_int);
    x_new = max(x_new-mult*dll,x_min);
    %%%%
    end
    t_iter(iter,i) = toc;
    fprintf('It.:%5i Vol.:%4.3f Eig1.:%7.3f Eig2.:%7.3f Eig3.:%7.3f mu:%5.3f L:%7.3f \n',...
        iter-1,(dll'*x),lambda(1),lambda(2),lambda(3),mu,L); 
    x = x_new;
end
[~] = draw_cs_freq(coord_x,irr,x,i,nx,ny);

% <--- 新增：计算并打印当前算法的总时间 --->
t_alg_total = toc(t_alg_start); 
fprintf('\n>>> 算法 %s 执行完毕！总耗时: %.2f [s] <<<\n\n', alg_list(i), t_alg_total);

end
%% Plot objective values
%x = x_new;
%end
%[~] = draw_cs_freq(coord_x,irr,x,i,nx,ny);
%end
%% Plot objective values
k = 0:maxiter;
figure(5)
for i = alg
    if i == 1
        plot(k,-data_obj(:,1),'-','linewidth',2,'Color',"#0072BD")
    elseif i == 2
        plot(k,-data_obj(:,2),'--','linewidth',2,'Color',"#D95319")
    elseif i == 3
        plot(k,-data_obj(:,3),'-.','linewidth',2,'Color',"#EDB120")
    elseif i == 4
        plot(k,-data_obj(:,4),':','linewidth',2,'Color',"#7E2F8E")
    end
    hold on
end
xlim([0,maxiter]);
xlabel('Iteration $k$','FontSize',15,'interpreter','latex');
ylabel('Objective value','FontSize',15,'interpreter','latex');
legend(alg_name,'interpreter','latex','Location','northeast')
hold off
set(gca,'FontSize',25)
ax = gca;
ax.FontSize = 22;
ax.TickLabelInterpreter = 'latex';
%% Plot objective values (log)
% opt_value = 2569.9479039758; % len_scale = 0.1, V0 = 1
% opt_value = 25.6938363532; % len_scale = 0.01, V0 = 1
% opt_value = 5109.0289540340; % len_scale = 1, V0 = 0.1, dmm0 = 10^(-2)
% opt_value = 51.2259438465; % len_scale = 1, V0 = 0.1, dmm0 = 1
opt_value = 51.4026958936; % aligned with main_bisection.m, nx=4, ny=4, Inf, bottom center
% opt_value =15.3984958306
% opt_value = 12.4998317912; % len_scale = 1, V0 = 0.1, dmm0 = 1, triple
% opt_value = 20.0239123041; % top center 19.9887729440
data_obj_diff = opt_value-data_obj;
k = 0:maxiter;
%% Plot obj
% figure(6)
% for i = alg
%     if i == 1
%         semilogy(k,data_obj_diff(:,1),'-','linewidth',2,'Color',"#0072BD")
%     elseif i == 2
%         semilogy(k,data_obj_diff(:,2),'--','linewidth',2,'Color',"#D95319")
%     elseif i == 3
%         semilogy(k,data_obj_diff(:,3),'-.','linewidth',2,'Color',"#EDB120")
%     elseif i == 4
%         semilogy(k,data_obj_diff(:,4),':','linewidth',2,'Color',"#7E2F8E")
%     end
%     hold on
% end
% xlim([0,maxiter]);
% ylim([0.0001,100]);
% xlabel('Iteration $k$','FontSize',15,'interpreter','latex');
% ylabel('$f(x^k)-f^*$','FontSize',15,'interpreter','latex');
% legend(alg_name,'interpreter','latex','Location','northeast')
% hold off
% set(gca,'FontSize',25)
% ax = gca;
% ax.FontSize = 22;
% ax.TickLabelInterpreter = 'latex';
%%%%
% fprintf('t_iter.: Alg1 = %7.6f Alg2 = %7.6f Alg3 = %7.6f Alg4 = %7.6f [s] \n',...
%     mean(t_iter(:,1)),mean(t_iter(:,2)),mean(t_iter(:,3)),mean(t_iter(:,4)));
% fprintf('Alg1 Eig1.:%10.7f Eig2.:%10.7f Eig3.:%10.7f \n',...
%     data_lambda(maxiter+1,1,1),data_lambda(maxiter+1,1,2),data_lambda(maxiter+1,1,3)); 
% fprintf('Alg2 Eig1.:%10.7f Eig2.:%10.7f Eig3.:%10.7f \n',...
%     data_lambda(maxiter+1,2,1),data_lambda(maxiter+1,2,2),data_lambda(maxiter+1,2,3)); 
% fprintf('Alg3 Eig1.:%10.7f Eig2.:%10.7f Eig3.:%10.7f \n',...
%     data_lambda(maxiter+1,3,1),data_lambda(maxiter+1,3,2),data_lambda(maxiter+1,3,3));
% fprintf('Alg4 Eig1.:%10.7f Eig2.:%10.7f Eig3.:%10.7f \n',...
%     data_lambda(maxiter+1,4,1),data_lambda(maxiter+1,4,2),data_lambda(maxiter+1,4,3)); 
% fprintf(' =========================================================== \n');
