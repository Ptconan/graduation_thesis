%%%%% 结构拓扑优化：三大等价问题统一数值验证脚本 (CVX 版)
%%%%% 用于在论文中证明 P_V, P_C, P_E 三者的全局最优拓扑图案完全等价
clear; clc;

% --- 1. 数据加载与基础参数 ---
nd = 46; 
n = 200; 
full_nd = 2*nd+1;
load('sqrtK.txt'); 
load('strmM.txt');
fl = zeros(nd, 1); 
fl(41) = 4;

% 预设的工程要求 (约束条件)
gamma_given = 1.0;     % 柔顺度要求
lambda_given = 5e6;    % 基频(特征值)要求

disp('======================================================');
disp('>> 开始求解 问题 1 (P_V): 最小化结构体积');
disp('约束条件: 柔顺度 <= 1.0, 基频 >= 5e6');
disp('======================================================');
tic;
cvx_begin sdp quiet
    variable x1(n, 1)
    minimize(sum(x1))
    subject to
        [gamma_given -fl' zeros(1,full_nd-nd-1); 
         -fl sqrtK*sparse(diag(x1))*sqrtK' zeros(nd,full_nd-nd-1); 
         zeros(full_nd-nd-1,nd+1) sqrtK*sparse(diag(x1))*sqrtK'-lambda_given*reshape(strmM*x1,nd,nd)] >= 0;
        x1 >= 0;
cvx_end
toc;
V_opt = sum(x1);
fprintf('【P_V 结果】最优体积 V* = %.4f\n\n', V_opt);


disp('======================================================');
disp('>> 开始求解 问题 2 (P_C): 最小化结构柔顺度');
fprintf('约束条件: 总体积 <= %.4f, 基频 >= 5e6\n', V_opt);
disp('======================================================');
tic;
cvx_begin sdp quiet
    variable x2(n, 1)
    variable c_comp
    minimize(c_comp)
    subject to
        [c_comp -fl' zeros(1,full_nd-nd-1); 
         -fl sqrtK*sparse(diag(x2))*sqrtK' zeros(nd,full_nd-nd-1); 
         zeros(full_nd-nd-1,nd+1) sqrtK*sparse(diag(x2))*sqrtK'-lambda_given*reshape(strmM*x2,nd,nd)] >= 0;
        sum(x2) <= V_opt;
        x2 >= 0;
cvx_end
toc;
fprintf('【P_C 结果】最优柔顺度 C* = %.4f (理论上应等于 %.4f)\n\n', c_comp, gamma_given);


disp('======================================================');
disp('>> 开始求解 问题 3 (P_E): 最大化最小特征值 (二分法搜索)');
fprintf('约束条件: 总体积 <= %.4f, 柔顺度 <= %.4f\n', V_opt, gamma_given);
disp('======================================================');
tic;
lambda_low = 1e5;  
lambda_high = 1e7; 
tol = 1e3;         

while (lambda_high - lambda_low) > tol
    lambda_mid = (lambda_low + lambda_high) / 2;
    
    cvx_begin sdp quiet
        variable x3_temp(n, 1)
        minimize(0) % 仅做可行性测试
        subject to
            [gamma_given -fl' zeros(1,full_nd-nd-1); 
             -fl sqrtK*sparse(diag(x3_temp))*sqrtK' zeros(nd,full_nd-nd-1); 
             zeros(full_nd-nd-1,nd+1) sqrtK*sparse(diag(x3_temp))*sqrtK'-lambda_mid*reshape(strmM*x3_temp,nd,nd)] >= 0;
            sum(x3_temp) <= V_opt;
            x3_temp >= 0;
    cvx_end
    
    if strcmp(cvx_status, 'Solved')
        lambda_low = lambda_mid;
        x3 = x3_temp; % 保存当前可行的拓扑解
    else
        lambda_high = lambda_mid;
    end
end
toc;
fprintf('【P_E 结果】最大可达基频 lambda* = %.2e (理论上应逼近 %.2e)\n\n', lambda_low, lambda_given);


disp('======================================================');
disp('>> 终极验证：拓扑结果 (向量 x) 数值等价性分析');
disp('======================================================');
% 计算相对误差 (Relative Error Norm)
err_12 = norm(x1 - x2) / norm(x1);
err_13 = norm(x1 - x3) / norm(x1);

fprintf('问题1 (P_V) 与 问题2 (P_C) 的拓扑结构相对误差: %.2e\n', err_12);
fprintf('问题1 (P_V) 与 问题3 (P_E) 的拓扑结构相对误差: %.2e\n', err_13);

if err_12 < 1e-3 && err_13 < 1e-3
    disp(' ');
    disp('【结论】证明成功！三个看似不同的数学问题，最终求得的各杆件截面积分布在数值上完全等价。');
else
    disp('【结论】误差偏大，可能是求解精度或求解器配置导致。');
end

% 可选：绘制对比图直观展示
figure('Name', 'Topology Variables Equivalence', 'Color', 'w');
plot(1:n, x1, 'bo', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(1:n, x2, 'r+', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(1:n, x3, 'k.', 'LineWidth', 1.5, 'MarkerSize', 10);
xlabel('杆件编号 (Bar Index)', 'FontSize', 12);
ylabel('截面积 (Cross-sectional Area)', 'FontSize', 12);
title('三个问题的拓扑设计变量重合度对比', 'FontSize', 14);
legend('问题 1 (极小化体积)', '问题 2 (极小化柔顺度)', '问题 3 (极大化基频)');
grid on;