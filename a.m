% 1. 设置与你论文一致的网格参数 (例如 4x4)
nx = 4; 
ny = 4;

% 2. 调用你的 member 函数生成节点和拓扑连接关系
[dll, matH, coord_x, ir, irr, ird, Idx] = member(nx, ny, Inf);
nm = size(matH, 2); % 获取总杆件数

% 3. 构造一个全为 1 的截面积向量（代表初始未优化的全连接基结构）
x_ground_structure = ones(nm, 1);

% 4. 调用你的画图函数，绘制出初始基结构
draw_cs_freq(coord_x, irr, x_ground_structure, 1, nx, ny);
title('初始密集基结构 (Ground Structure)');