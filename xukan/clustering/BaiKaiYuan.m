%白凯远的基于双向伪距测量的分布式时间同步算法性能仿真

%% 无人机蜂群运动模型
% 特点：5km级战场范围、500m通信半径、稀疏拓扑、高速机动
%history.adj_mat{t} ,history.dist_mat{t} 分别代表连接矩阵和距离矩阵
num_drones = 50;
steps = 2000;           % 维持长航时
dt = 0.1;

% --- 视觉参数 ---
plot_stride = 3;
view_radius = 1000;     % <--- 视野半径扩大到1km，以容纳庞大的稀疏蜂群

% --- 关键：空间参数 (按比例扩大至 500m 通信级) ---
comm_range = 500;       % <--- 通信半径 (目标值)
r_separation = 250;     % <--- 排斥半径 (保持 0.6*comm_range 比例，维持稀疏度)
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
% %% <--- 修改 1：初始化坐标记录矩阵 --->
% 维度：[时间步, 节点数, 3(x,y,z)]
history.pos = zeros(steps, num_drones, 3);

%% 3. 仿真循环
fprintf('启动超广域仿真 (通信半径: %dm)...\n', comm_range);



last_adj = zeros(num_drones);

for t = 1:steps
    % --- A. 拓扑计算 ---
    dist_mat = pdist2(pos, pos);
    adj_mat = (dist_mat < comm_range) & (dist_mat > 0);

    % 记录
    history.adj_mat{t} = double(adj_mat);
    history.dist_mat{t} = dist_mat;
    % %% <--- 修改 3：记录当前坐标 --->
    history.pos(t, :, :) = pos;

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

end
fprintf('仿真完成。\n');


%%利用MTS进行粗同步，同时进行簇头节点选举


t_period=1;%发送消息伪周期
time_step = 0.1; % 时间步长（s）
simulation_k=185;%仿真轮数
%不同轮次调整后的频偏和相偏
x=zeros(num_drones,simulation_k);
y=zeros(num_drones,simulation_k);
z=zeros(num_drones,simulation_k);
delta_max=zeros(simulation_k,1);%每轮迭代后全网最大钟差
t=zeros(simulation_k,1);%每轮迭代后用于测评的真实时间，此处选为最后一个发送消息的节点每次发送完的时间
skew_final=zeros(1,simulation_k);%平均后的每一轮频偏差值。
offset_final=zeros(1,simulation_k);
clock_final=zeros(1,simulation_k);
phy_skew=80e-6;%物理时钟频偏最大值
phy_offset=1;%物理时钟相位偏移最大值


%初始化物理时钟频偏
alpha=1-phy_skew+(2*phy_skew)*rand(num_drones,1);
%初始化物理时钟相偏
beta=-phy_offset+2*phy_offset*rand(num_drones,1);
%初始化逻辑时钟频率调整参数
l=ones(num_drones,1);
%初始化逻辑时钟相位调整参数
h=zeros(num_drones,1);


t_local_total=zeros(num_drones,simulation_k);
t_global_total=zeros(num_drones,simulation_k);
t_global=(t_period*0.9)*sort(rand(num_drones,1));%初始发送全球时间，升序,预留0.1倍周期作为保护jiange
t_local=alpha.*t_global+beta;%初始发送节点时间
t_local_total(:,1)=t_local;
t_global_total(:,1)=t_global;

%所有轮次中各节点发送消息的本地物理时间
for m=2:simulation_k
    t_local_total(:,m)=t_local+t_period*(m-1);
end

%所有轮次中各节点发送消息的全局时间
for m=2:simulation_k
    t_global_total(:,m)=(t_local_total(:,m)-beta)./alpha;
end

%初始化每个节点本地的邻居调整信息列表(消息轮，约定发送物理时间，自身接收物理，邻居接收物理，上一时刻的ηij,收到同一节点总轮次)，
%前n行表示首次接收到某节点的消息数据
%后n行表示最新接收到某节点的消息数据
neb_list=cell(num_drones,1);
for i=1:num_drones
    neb_list{i}=zeros(num_drones*2,5);
end

for j=1:simulation_k
    for i=1:num_drones
        [~,d]=sort(t_global_total(:,j));

        A=history.adj_mat{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的连接矩阵
        B1=history.dist_mat{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的距离矩阵

        %第一轮不更新，从第二轮开始
        if j>=2
            %(消息轮，约定发送物理时间，自身接收物理，邻居接收物理，上一时刻的alphaji,收到同一节点总轮次)，
            self_list=neb_list{d(i)};%自身的调整矩阵
            neb_count_l=0;%邻居数量
            neb_count_h=0;
            l_change=0;
            h_change=0;
            for ii=1:num_drones
                %遍历矩阵的下n行
                %如果有最新消息，则利用两轮消息计算相对频偏和时钟差值，更新上一时刻相对频偏值
                if self_list(ii+num_drones,1)>=j-1
                    rou_ji_k=self_list(ii,3)-self_list(ii,2);
                    rou_ij_k=self_list(ii,4)-self_list(ii,2);
                    rou_ji_k1=self_list(ii+num_drones,3)-self_list(ii+num_drones,2);
                    rou_ij_k1=self_list(ii+num_drones,4)-self_list(ii+num_drones,2);
                    S_ij=(rou_ji_k1-rou_ij_k1)/2-(rou_ji_k-rou_ij_k)/2;
                    alpha_ji=S_ij/(self_list(ii+num_drones,2)-self_list(ii,2))+...
                        sqrt(1+(S_ij/(self_list(ii+num_drones,2)-self_list(ii,2)))^2);
                    if self_list(ii,5)~=0
                        alpha_ji=self_list(ii,5)/2+alpha_ji/2;
                        self_list(ii,5)=alpha_ji;
                    else
                        self_list(ii,5)=alpha_ji;
                    end
                    neb_list{d(i)}=self_list;
                    neb_count_l=neb_count_l+1;
                    neb_count_h=neb_count_h+1;
                    l_change=l_change+(1/alpha_ji)*l(ii,1)-l(d(i));
                    h_change=h_change+l(ii,1)*self_list(ii+num_drones,2)+h(ii,1)-...
                        l(d(i),1)*(self_list(ii+num_drones,3)-B1(ii,d(i))/3e8)-h(d(i),1);
                    flood_range=(rou_ij_k1+rou_ji_k1)/2;
                    flood_range_real=B1(ii,d(i))/3e8;
                    %alpha_ji_real=alpha(d(i),1)/alpha(ii,1);

                    %如果没有最新消息，检查上n行是否有最新消息，如果有，利用一轮消息计算时钟差值
                elseif self_list(ii,1)>=j-1
                    rou_ji_k=self_list(ii,3)-self_list(ii,2);
                    rou_ij_k=self_list(ii,4)-self_list(ii,2);


                    %flood_range(ii,d(i))=(rou_ij_k+rou_ji_k)/2;
                    h_change=h_change+l(ii,1)*self_list(ii,2)+h(ii,1)-...
                        l(d(i),1)*(self_list(ii,3)-B1(ii,d(i))/3e8)-h(d(i),1);
                    neb_count_h=neb_count_h+1;
                end
            end
            %频率补偿参数调整
            l(d(i),1)=l(d(i),1)+l_change/(neb_count_l+1);
            %相位补偿参数调整
            h(d(i),1)=h(d(i),1)+h_change/(neb_count_h+1);
        end


        neb_count=0;%用于确定和每个邻居同时发送的约定时间
        %检查d(i)节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
        for k=1:num_drones
            if A(d(i),k)==1

                F=neb_list{k};%邻居信息矩阵

                %首次接收到的消息（不管是第几轮），存放在前n行，n表示无人机总数
                if F(d(i),1)==0
                    F(d(i),1)=j;
                    %填充邻居列表中首次接收到该节点信息的数据
                    F(d(i),2)=t_local_total(d(i),j)+0.01*neb_count+5e-9*randn;%约定的同时发送时间

                    %邻居接收到的时间
                    F(d(i),3)=(((t_local_total(d(i),j)+0.01*neb_count-...
                        beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;

                    %自己接收到的时间+100*rand-50
                    F(d(i),4)=(((t_local_total(d(i),j)+0.01*neb_count-...
                        beta(k,1))/(alpha(k,1)))+abs(B1(d(i),k)+(70*rand-35)*(beta(k,1)-beta(d(i),1)))/3e8)*alpha(d(i),1)+beta(d(i),1)+5e-9*randn;

                    %flood_range(d(i),k)=(F(d(i),3)-F(d(i),2)+F(d(i),4)-F(d(i),2))/2;
                    neb_list{k}=F;
                    neb_count=neb_count+1;

                    %第二次接受到的消息（不管是那一轮），
                    %频率消息存放在第n+1到2n行，用于计算的数据是首次数据和最新数据
                elseif F(d(i)+num_drones,1)==0
                    F(d(i)+num_drones,1)=j;

                    %填充邻居列表中首次接收到该节点信息的数据
                    F(d(i)+num_drones,2)=t_local_total(d(i),j)+0.01*neb_count+5e-9*randn;%约定的同时发送时间

                    %邻居接收到的时间
                    F(d(i)+num_drones,3)=(((t_local_total(d(i),j)+0.01*neb_count-...
                        beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;

                    %自己接收到的时间
                    F(d(i)+num_drones,4)=(((t_local_total(d(i),j)+0.01*neb_count-...
                        beta(k,1))/(alpha(k,1)))+abs(B1(d(i),k)+(70*rand-35)*(beta(k,1)-beta(d(i),1)))/3e8)*alpha(d(i),1)+beta(d(i),1)+5e-9*randn;
                    neb_list{k}=F;
                    neb_count=neb_count+1;


                    %第三次或以上接受到的消息（不管是那一轮），
                    %n+1-2n数据放回前n行
                    %频率消息存放在第n+1到2n行，用于计算的数据是上次数据和最新数据
                else
                    %数据网上覆盖
                    F(d(i),1:4)=F(d(i)+num_drones,1:4);
                    %更新信息矩阵
                    F(d(i)+num_drones,1)=j;

                    %填充邻居列表中首次接收到该节点信息的数据
                    F(d(i)+num_drones,2)=t_local_total(d(i),j)+0.01*neb_count+5e-9*randn;%约定的同时发送时间

                    %邻居接收到的时间
                    F(d(i)+num_drones,3)=(((t_local_total(d(i),j)+0.01*neb_count-...
                        beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;

                    %自己接收到的时间
                    F(d(i)+num_drones,4)=(((t_local_total(d(i),j)+0.01*neb_count-...
                        beta(k,1))/(alpha(k,1)))+abs(B1(d(i),k)+(70*rand-35)*(beta(k,1)-beta(d(i),1)))/3e8)*alpha(d(i),1)+beta(d(i),1)+5e-9*randn;
                    neb_list{k}=F;
                    neb_count=neb_count+1;
                end
            end
        end

    end
    %最终调整的逻辑时钟频偏
    x(:,j)=l.*alpha;
    %最终调整的逻辑时钟相偏
    y(:,j)=l.*beta+h;
    z(:,j)=h;
end
%选取每次评估同步效果的时间（每轮结束评估一次）
for i=1:simulation_k
    t(i,1)=(t_local_total(num_drones,i)-beta(num_drones))/alpha(num_drones);
end
%每轮结束后的时钟偏差
for k=1:simulation_k
    delta_max(k,1)=max(x(:,k)*t(k,1)+y(:,k))-min(x(:,k)*t(k,1)+y(:,k));
end


x1=zeros(simulation_k,1);
y1=zeros(simulation_k,1);
%每轮结束后的频率偏差
for c=1:simulation_k
    x1(c,1)=max(x(:,c))-min(x(:,c));
end
%每轮结束后的相位偏差
for d=1:simulation_k
    y1(d,1)=max(y(:,d))-min(y(:,d));
end
skew_final=x1';
offset_final=y1';
clock_final=delta_max';




format long;
% semilogy(1:5:simulation_k,skew_final(1:5:end),'-*','LineWidth',1.5);
% hold on;
% semilogy(1:5:simulation_k,offset_final(1:5:end),'-*','LineWidth',1.5);
%
% hold on;
semilogy(1:5:simulation_k,clock_final(1:5:end),'-*','LineWidth',1.5);
hold on;
xlabel("同步轮数");
ylabel("全网最大频率偏移");
title("MTS粗同步频率偏差")

% semilogy(1:simulation_k,clock_var_final(1:end),'-*','LineWidth',1);
grid on;

aaaaa=1-1/(1+8e-7);