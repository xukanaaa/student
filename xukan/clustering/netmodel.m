%% 超广域战术蜂群仿真 (Tactical Scale Swarm)
% 特点：5km级战场范围、500m通信半径、稀疏拓扑、高速机动
% 视觉：白色背景 + 浅色连线 + 广角镜头跟随

clear; clc; close all;

%% 1. 参数全比例放大 (Scale Up)
num_drones = 50;        
steps = 2000;           % 维持长航时
dt = 0.1;               

% --- 视觉参数 ---
plot_stride = 3;        
view_radius = 1000;     % <--- 视野半径扩大到1km，以容纳庞大的稀疏蜂群

% --- 关键：空间参数 (按比例扩大至 500m 通信级) ---
comm_range = 500;       % <--- 通信半径 (目标值)
r_separation = 300;     % <--- 排斥半径 (保持 0.6*comm_range 比例，维持稀疏度)
max_speed = 50;         % <--- 最大速度 (提高以适应大地图)

% --- 动力学权重 (针对大尺度距离进行补偿) ---
w_target     = 0.8;     
w_separation = 15.0;    % <--- 大幅增加权重，补偿距离平方衰减导致的斥力微弱
w_alignment  = 0.4;     
w_cohesion   = 0.1;     % 弱凝聚
w_noise      = 10.0;    % <--- 强随机扰动，在大尺度下制造肉眼可见的拓扑变化

%% 2. 初始化 (战术级分布)
% 初始生成范围扩大到 1500x1500
pos = rand(num_drones, 3) .* [1500, 1500, 400]; 
vel = (rand(num_drones, 3) - 0.5) * max_speed;

% 战术巡航路径 (坐标范围扩展到 5000m+)
waypoints = [
    1000, 1000, 400;
    3000, 1500, 600;   % 远距离爬升
    4500, 3000, 500;   % 深入腹地
    3500, 4500, 800;   % 高空盘旋侦察
    1500, 4000, 600;   % 侧向机动
    500,  2000, 300;   % 返航段
    1000, 1000, 400    % 闭环
];
current_wp_idx = 1;

% 数据存储
history.adj_mat = cell(steps, 1);
history.dist_mat = cell(steps, 1);

%% 3. 仿真循环
fprintf('启动超广域仿真 (通信半径: %dm)...\n', comm_range);

figure('Name', 'Tactical Swarm Simulation (500m Comm)', 'Color', 'w', 'Position', [50, 50, 1200, 800]);
h_ax = gca;
axis equal; grid on; hold on;
view(3);
xlabel('East (m)'); ylabel('North (m)'); zlabel('Altitude (m)');

% --- 绘图风格 ---
% 节点：深蓝，稍微画大一点点以便在大视野下看清
h_nodes = scatter3(0,0,0, 30, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'none'); 

% 连线：浅灰色 [0.75 0.75 0.75] + 透明度
h_links = plot3([0,0],[0,0],[0,0], 'Color', [0.75 0.75 0.75 0.5], 'LineWidth', 0.8); 

% 路径：红色虚线
plot3(waypoints(:,1), waypoints(:,2), waypoints(:,3), 'r--', 'LineWidth', 1.2); 
title_h = title('Initializing...');

last_adj = zeros(num_drones);

for t = 1:steps
    % --- A. 拓扑计算 ---
    dist_mat = pdist2(pos, pos);
    adj_mat = (dist_mat < comm_range) & (dist_mat > 0);
    
    % 记录
    history.adj_mat{t} = double(adj_mat);
    history.dist_mat{t} = dist_mat;
    
    % 统计
    topo_changes = sum(sum(xor(adj_mat, last_adj))) / 2;
    last_adj = adj_mat;
    avg_deg = mean(sum(adj_mat));
    
    % --- B. 物理受力 (Physics Scale-Up) ---
    target = waypoints(current_wp_idx, :);
    
    % 切换路点判定距离放大
    if norm(mean(pos) - target) < 400 
        current_wp_idx = mod(current_wp_idx, size(waypoints,1)) + 1;
    end
    
    sep_force = zeros(num_drones, 3);
    coh_force = zeros(num_drones, 3);
    align_force = zeros(num_drones, 3);
    
    for i = 1:num_drones
        d = dist_mat(i, :);
        % 感知范围
        neighbors = find(d < comm_range * 1.2 & d > 0); 
        
        if ~isempty(neighbors)
            % 分离力
            close_mask = d(neighbors) < r_separation;
            for j = neighbors(close_mask)
                diff = pos(i,:) - pos(j,:);
                dist_val = norm(diff);
                % 重要：距离平方项在大尺度下数值巨大(300^2=90000)，导致力微乎其微
                % 修正：在分母中除以一个比例因子(例如1000)或者直接依靠增大的权重
                sep_force(i,:) = sep_force(i,:) + (diff / (dist_val^2 + 10.0));
            end
            % 凝聚与对齐
            coh_force(i,:) = mean(pos(neighbors, :), 1) - pos(i,:);
            align_force(i,:) = mean(vel(neighbors, :), 1) - vel(i,:);
        else
            % 孤立回归
            coh_force(i,:) = (mean(pos) - pos(i,:)) * 5; 
        end
    end
    
    % 归一化
    norm_v = @(v) v ./ (sqrt(sum(v.^2, 2)) + 1e-6) * max_speed;
    noise = (rand(num_drones, 3) - 0.5) * max_speed;
    
    % 目标力
    t_vec = target - pos;
    t_force = (t_vec ./ (sqrt(sum(t_vec.^2,2)) + 1e-6)) * max_speed;
    
    % 合成力 (注意这里的乘数调整)
    total_force = w_target * t_force + ...
                  w_separation * sep_force * 100 + ... % <--- 乘数加大到100，抵消距离平方带来的数值衰减
                  w_alignment * norm_v(align_force) + ...
                  w_cohesion * norm_v(coh_force) + ...
                  w_noise * noise;
              
    vel = vel + total_force * dt;
    % 限速
    s = sqrt(sum(vel.^2, 2));
    idx = s > max_speed;
    vel(idx, :) = (vel(idx, :) ./ s(idx)) * max_speed;
    
    pos = pos + vel * dt;
    pos(:,3) = max(pos(:,3), 100); % 最低飞行高度 100m
    
    % --- C. 可视化更新 ---
    if mod(t, plot_stride) == 0
        % 更新节点
        set(h_nodes, 'XData', pos(:,1), 'YData', pos(:,2), 'ZData', pos(:,3));
        
        % 更新连线
        [r, c] = find(triu(adj_mat));
        if ~isempty(r)
            lx = [pos(r,1) pos(c,1) nan(size(r))]';
            ly = [pos(r,2) pos(c,2) nan(size(r))]';
            lz = [pos(r,3) pos(c,3) nan(size(r))]';
            delete(h_links);
            % 浅灰色连线
            h_links = line(h_ax, lx(:), ly(:), lz(:), 'Color', [0.75 0.75 0.75 0.5], 'LineWidth', 0.8);
        end
        
        % 镜头跟随 (view_radius 已扩大到 1000m)
        center = mean(pos, 1);
        xlim([center(1)-view_radius, center(1)+view_radius]);
        ylim([center(2)-view_radius, center(2)+view_radius]);
        zlim([max(0, center(3)-view_radius/2), center(3)+view_radius/2]);
        
        % 标题
        set(title_h, 'String', sprintf('Step: %d | Avg Neighbors: %.1f | Comm: %dm | Changes: %d', ...
            t, avg_deg, comm_range, topo_changes));
        
        drawnow;
    end
end
fprintf('仿真完成。\n');