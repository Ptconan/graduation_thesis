%%%%% 论文级配图生成：CVX 与 BM-ALM 在不同规模下的性能对比
%%%%% 注：本代码使用模拟的趋势数据生成模板，撰写论文时请将横纵坐标替换为你实际测出的真实数据！
clear; clc; close all;

% ==================== 1. 准备实验数据 ====================
% 假设我们测试了 8 组不同规模的网格 (变量数 N)
N_scale = [200, 500, 1000, 2000, 3000, 4000, 5000, 8000, 10000];

% 模拟 CVX(内点法) 的内存消耗 (O(N^2) 爆炸式增长，单位 MB)
mem_cvx = 0.001 * N_scale.^2 + 50; 
% 假设电脑内存上限为 16GB (约 16000MB)，N>=4000 时崩溃，数据变为 NaN
mem_cvx(N_scale > 4000) = NaN; 

% 模拟 BM-ALM(一阶法) 的内存消耗 (O(N) 线性增长，单位 MB)
mem_bm = 0.08 * N_scale + 100; 

% 模拟时间消耗 (单位: 秒)
time_cvx = 1e-7 * N_scale.^3 + 0.5; % 内点法解海森矩阵时间 O(N^3)
time_cvx(N_scale > 4000) = NaN;     % 崩溃后无时间数据
time_bm = 0.005 * N_scale.^1.5 + 5; % BM-ALM 时间稳步增长

% ==================== 2. 绘制内存消耗对比图 ====================
fig1 = figure('Name', 'Memory Comparison', 'Position', [100, 100, 600, 450], 'Color', 'w');

% 绘制曲线 (使用适合打印的粗细和标记)
plot(N_scale, mem_cvx, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980]); hold on;
plot(N_scale, mem_bm, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410]);

% 绘制 16GB 内存极限辅助线
yline(16000, 'k--', 'LineWidth', 1.5);
text(5000, 16500, 'RAM Limit (e.g., 16GB Out of Memory)', 'FontSize', 11, 'FontName', 'Times New Roman', 'Color', 'k');

% 图表美化设置
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'LineWidth', 1.2);
xlabel('Number of Variables (Scale $N$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Memory Consumption (MB)', 'Interpreter', 'latex', 'FontSize', 14);
legend('CVX (Interior Point Method)', 'BM-ALM (First-order Method)', 'Location', 'northwest', 'FontSize', 11);
grid on; box on;

% 导出为高质量 PDF 矢量图
exportgraphics(fig1, 'Fig_Memory_Comparison.pdf', 'ContentType', 'vector');
disp('成功导出内存对比矢量图：Fig_Memory_Comparison.pdf');


% ==================== 3. 绘制计算耗时对比图 ====================
fig2 = figure('Name', 'Time Comparison', 'Position', [750, 100, 600, 450], 'Color', 'w');

plot(N_scale, time_cvx, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980]); hold on;
plot(N_scale, time_bm, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410]);

% 图表美化设置
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'LineWidth', 1.2);
xlabel('Number of Variables (Scale $N$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Computational Time (Seconds)', 'Interpreter', 'latex', 'FontSize', 14);
legend('CVX (Interior Point Method)', 'BM-ALM (First-order Method)', 'Location', 'northwest', 'FontSize', 11);
grid on; box on;

% 导出为高质量 PDF 矢量图
exportgraphics(fig2, 'Fig_Time_Comparison.pdf', 'ContentType', 'vector');
disp('成功导出耗时对比矢量图：Fig_Time_Comparison.pdf');