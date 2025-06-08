num_drones = 100; % 无人机数量
space_size = 3000; % 空间大小（长、宽、高均为1km）
min_speed = 80; % 最小速度（m/s）
max_speed = 100; % 最大速度（m/s）
comm_range = 00; % 通信范围（m）
time_step = 1; % 时间步长（s）
simulation_time = 100; % 模拟时间（s）
speed_change=zeros(num_drones,3);%无人机初始化速度后相邻时间（1秒）速度的变化
% 初始化无人机位置
drone_positions = space_size * rand(num_drones, 3);
% 初始化无人机速度
drone_speeds = (max_speed - min_speed) * rand(num_drones, 3) + min_speed;
% 确保速度向量的模在4-8m/s范围内（调整方向但不改变模）
drone_speeds = drone_speeds ./ sqrt(sum(drone_speeds.^2, 2)) .* (min_speed + (max_speed - min_speed) * rand(num_drones, 1));

% 初始化无人机通信邻接矩阵
comm_matrix = zeros(num_drones);

%初始化每个离散时间点的距离矩阵
for i=1:(simulation_time/time_step)
    D{i}=zeros(num_drones,num_drones);
end
%初始化每个离散时间点的连接矩阵
for i=1:(simulation_time/time_step)
    I{i}=zeros(num_drones,num_drones);
end

% 模拟无人机运动
for k = 1:simulation_time/time_step
    speed_change=10*randn(num_drones,3);
    drone_speeds=drone_speeds+speed_change;%更新速度的扰动
    % 更新无人机位置
    drone_positions = drone_positions + drone_speeds * time_step;
    
    % 确保无人机不会超出空间边界
    for i=1:num_drones
        if drone_positions(i,1)<0
            drone_positions(i,1)=-drone_positions(i,1);
            drone_speeds(i,1)=-drone_speeds(i,1);
        end
        if drone_positions(i,1)>space_size
            drone_positions(i,1)=2*space_size-drone_positions(i,1);
            drone_speeds(i,1)=-drone_speeds(i,1);
        end
        if drone_positions(i,2)<0
            drone_positions(i,2)=-drone_positions(i,2);
            drone_speeds(i,2)=-drone_speeds(i,2);
        end
        if drone_positions(i,2)>space_size
            drone_positions(i,2)=2*space_size-drone_positions(i,2);
            drone_speeds(i,2)=-drone_speeds(i,2);
        end
        if drone_positions(i,3)<0
            drone_positions(i,3)=-drone_positions(i,3);
            drone_speeds(i,3)=-drone_speeds(i,3);
        end
        if drone_positions(i,3)>space_size
            drone_positions(i,3)=2*space_size-drone_positions(i,3);
            drone_speeds(i,3)=-drone_speeds(i,3);
        end
    end
    
    % 计算无人机之间的距离
    distances = pdist2(drone_positions, drone_positions);
    D{k}=distances;
    
    % 更新通信邻接矩阵
    comm_matrix = (distances <= comm_range);
    comm_matrix = comm_matrix-eye(num_drones,num_drones);
    I{k}=comm_matrix;
    
    %可视化无人机位置和通信关系（可选）
        figure;
        scatter3(drone_positions(:,1), drone_positions(:,2), drone_positions(:,3), 'filled');
        hold on;
        for i = 1:num_drones
            for j = i+1:num_drones
                if comm_matrix(i,j)
                    line([drone_positions(i,1), drone_positions(j,1)], ...
                         [drone_positions(i,2), drone_positions(j,2)], ...
                         [drone_positions(i,3), drone_positions(j,3)], ...
                         'Color', 'r', 'LineWidth', 1);
                end
            end
        end
        grid on;
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title(['Simulation Time: ', num2str(k), ' s']);
        drawnow;
    d=D{k};
    i=I{k};
    disp(d(1,2));
        disp(i(1,2));
    %暂停一段时间以便观察（可选）
        pause(5);
end


