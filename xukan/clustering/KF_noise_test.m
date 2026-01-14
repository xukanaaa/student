
%测试卡尔曼滤波对频偏的跟踪效果
%模型简化为两个节点间的同步参数估计

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


%===========================================================================
%%利用MTS进行粗同步，同时进行簇头节点选举
%1-5:5，6-15：10,16-30:15,31-50：20

t_period=1;%发送消息伪周期
time_step = 0.1; % 时间步长（s）
simulation_k=80;%仿真轮数
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

%初始化每个节点本地的邻居频率调整信息列表(消息轮，(当前邻居物理，当前自身物理)，当前邻居逻辑时钟频率调整l,
%flag位，表示当前同步状态，连续两轮自己的不变，且邻居的也不变才为1，变化为0，一轮为0.5，二轮为1
%前n行表示首次接收到某节点的消息数据
%后n行表示最新接收到某节点的消息数据
neb_list_fre=cell(num_drones,1);
for i=1:num_drones
    neb_list_fre{i}=zeros(num_drones*2,5);
end

%初始化每个节点本地的相位邻居信息列表(消息轮，(当前邻居物理，当前自身物理)，邻居逻辑时钟频率调整l，邻居逻辑时钟相位调整h
%,)
%前n行表示首次接收到某节点的消息数据
neb_list_pha=cell(num_drones,1);
for i=1:num_drones
    neb_list_pha{i}=zeros(num_drones,5);
end

%初始化每个节点收到的头结点得分汇总
neb_list_score=cell(num_drones,1);
for i=1:num_drones
    neb_list_score{i}=zeros(num_drones,1);
end

round_cluster=0;

for j=1:simulation_k
    for i=1:num_drones
        [~,d]=sort(t_global_total(:,j));

        %第一轮不更新，从第二轮开始
        if j>=2
            M=mts_l_h(neb_list_fre{d(i)},neb_list_pha{d(i)},j,num_drones,l(d(i),1),h(d(i),1),d(i));%更新逻辑时钟频偏和相偏参数
            %返回值包括频偏补偿，相偏补偿，当前同步状态（0,0.5,1）
            l(d(i),1)=M(1,1);
            h(d(i),1)=M(1,2);
            state_self=neb_list_fre{d(i)};
            state_self(d(i),5)=M(1,3);
            neb_list_fre{d(i)}=state_self;

        end
        A=history.adj_mat{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的连接矩阵
        B1=history.dist_mat{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的距离矩阵

        %更新完后对当前节点局部区域同步状态进行判断，如果自己连续两轮未改变且当前轮次收到的邻居也是如此，状态为1，否则为0
        F1=neb_list_fre{d(i)};%广播节点的自身的消息矩阵
        syn_state=1;
        for ii=num_drones+1:2*num_drones
            if F1(ii,1)>=j-1
                if F1(ii,5)<1
                    syn_state=0;
                end
            end
        end
        if F1(d(i),5)<1
            syn_state=0;
        end

        %计算当前节点score，当且仅当syn_state为1，且自己的S矩阵对应的自己的值为0时才计算
        S_self=neb_list_score{d(i)};
        if S_self(d(i),1)==0&&syn_state==1
            score_self=get_score(d(i),A,num_drones,j);
            S_self(d(i),1)=score_self;
        end

        %计算簇头选举出来时的轮次
        if all(S_self(:))
            if S_self(d(i),1)==min(S_self)&&round_cluster==0
                round_cluster=j;%簇头选举出来的轮次
                max_node=d(i);
                A_cluster=A;%连接矩阵
                B_cluster=B1;%距离矩阵
                S_self_temp=S_self;
                fre_initial=l.*alpha;
                pha_initial=l.*beta+h;
            end
        end

        %检查d(i)节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
        for k=1:num_drones
            if A(d(i),k)==1

                F=neb_list_fre{k};%频率矩阵
                P=neb_list_pha{k};%相位矩阵
                S=neb_list_score{k};%邻居的分数矩阵

                %首次接收到的消息（不管是第几轮），存放在前n行，n表示无人机总数
                if F(d(i),1)==0
                    F(d(i),1)=j;
                    P(d(i),1)=j;
                    %填充邻居列表中首次接收到该节点信息的数据(频率和相位矩阵都有)
                    F(d(i),2)=t_local_total(d(i),j)+5e-9*randn;
                    P(d(i),2)=F(d(i),2);

                    F(d(i),3)=(((t_local_total(d(i),j)-...
                        beta(d(i),1))/(alpha(d(i),1))))*alpha(k,1)+beta(k,1)+2*B1(d(i),k)/3e8+5e-9*randn;
                    P(d(i),3)=F(d(i),3);

                    F(d(i),4)=l(d(i),1);
                    P(d(i),4)=l(d(i),1);
                    P(d(i),5)=h(d(i),1);

                    %更新合并自己的分数矩阵和邻居的分数矩阵
                    S_merge=merge(S,S_self,num_drones);
                    S=S_merge;

                    neb_list_fre{k}=F;
                    neb_list_pha{k}=P;
                    neb_list_score{k}=S;
                end

                %非首次最新接受到的消息（不管是那一轮），
                %频率消息存放在第n+1到2n行，用于计算的数据是首次数据和最新数据
                %相位消息更新
                if F(d(i),1)~=0
                    F(d(i)+num_drones,1)=j;
                    P(d(i),1)=j;

                    F(d(i)+num_drones,2)=t_local_total(d(i),j)+5e-9*randn;


                    F(d(i)+num_drones,3)=(((t_local_total(d(i),j)-...
                        beta(d(i),1))/(alpha(d(i),1))))*alpha(k,1)+beta(k,1)+2*B1(d(i),k)/3e8+5e-9*randn;



                    %更新相位矩阵的值
                    P(d(i),2)=F(d(i)+num_drones,2);
                    P(d(i),3)=F(d(i)+num_drones,3);



                    P(d(i),4)=l(d(i),1);
                    P(d(i),5)=h(d(i),1);
                    F(d(i)+num_drones,4)=l(d(i),1);

                    F_self=neb_list_fre{d(i)};
                    F(d(i)+num_drones,5)=F_self(d(i),5);

                    %更新合并自己的分数矩阵和邻居的分数矩阵
                    S_merge=merge(S,S_self,num_drones);
                    S=S_merge;

                    neb_list_fre{k}=F;
                    neb_list_pha{k}=P;
                    neb_list_score{k}=S;
                end
            end
        end
    end
    %最终调整的逻辑时钟频偏
    x_mts(:,j)=l.*alpha;
    %最终调整的逻辑时钟相偏
    y_mts(:,j)=l.*beta+h;
end
%选取每次评估同步效果的时间（每轮结束评估一次）
for i=1:simulation_k
    t(i,1)=(t_local_total(num_drones,i)-beta(num_drones))/alpha(num_drones);
end
%每轮结束后的时钟偏差
for k=1:simulation_k
    delta_max_mts(k,1)=max(x_mts(:,k)*t(k,1)+y_mts(:,k))-min(x_mts(:,k)*t(k,1)+y_mts(:,k));
end

x1_mts=zeros(simulation_k,1);
y1_mts=zeros(simulation_k,1);
%每轮结束后的频率偏差
for c=1:simulation_k
    x1_mts(c,1)=max(x_mts(:,c))-min(x_mts(:,c));
end
%每轮结束后的相位偏差
for d=1:simulation_k
    y1_mts(d,1)=max(y_mts(:,d))-min(y_mts(:,d));
end
%一共进行5次蒙特卡洛仿真，r=5，每一行表示各轮次同步后的全网最大频率差
skew_mts=x1_mts';
offset_mts=y1_mts';
clock_mts=delta_max_mts';






%卡尔曼测试部分，这里只选取1号节点和二号节点进行同步
%中间发生频率的跳变。
%比较卡尔曼滤波和最小二乘线性回归，低通滤波对突变值的跟踪
%=====================================================
flood_alpha=fre_initial;%泛洪同步频率初始值
flood_beta=pha_initial;%泛洪同步相位初始值
t_master=flood_alpha(1,1)+flood_beta(1,1);%主节点初次泛洪时的本地物理时间
l_KF=ones(200,1);%表示每一轮卡尔曼的补偿值
l_LR=ones(200,1);%表示每一轮线性回归的补偿值
l_LP=ones(200,1);%表示每一轮低通滤波的补偿值


%% 卡尔曼滤波

%初始化2节点收到的信息
%每个矩阵元素代表：泛洪轮次，节点上次发送时间，上次自己接收时间，上次对方接收时间，卡尔曼滤波alpha，卡尔曼滤波beta
%最后四列表示协方差矩阵11,12,21,22位置的值
neb_list_flood_KF=zeros(num_drones,10);
flood_l=ones(1,1);%泛洪频率补偿
flood_h=zeros(1,1);%泛洪相位补偿

%主节点以1s为周期与邻居进行同步
for round_KF=1:200



    %发送消息时刻的距离矩阵
    flood_B1=history.dist_mat{ceil((t_master-flood_beta(1,1))/(flood_alpha(1,1)*time_step))};
    t_send_i_j=t_master+5e-9*randn;%二者约定的发送时间
    %邻居接收时间计算
    t_receive_j=(((t_master-...
        flood_beta(1,1))/(flood_alpha(1,1)))+flood_B1(1,2)/3e8)*flood_alpha(2,1)+...
        flood_beta(2,1)+5e-9*randn;
    %自身接收时间计算
    t_receive_i=(((t_master-...
        flood_beta(2,1))/(flood_alpha(2,1)))+flood_B1(1,2)/3e8)*flood_alpha(1,1)+...
        flood_beta(1,1)+5e-9*randn;

    % 每个矩阵元素代表：泛洪轮次，节点上次发送时间，上次自己接收时间，上次对方接收时间，
    % 卡尔曼滤波alpha，卡尔曼滤波beta
    flood_neb_history_KF=neb_list_flood_KF;%邻居节点的历史信息矩阵

    %邻居节点第一次接收到该节点信息，只进行相位同步
    if flood_neb_history_KF(1,1)==0
        rou_ij=t_receive_j-t_send_i_j;
        rou_ji=t_receive_i-t_send_i_j;
        flood_delta_t=(rou_ji-rou_ij)/2;%j比i慢多少
        flood_range=(rou_ji+rou_ij)/2;%二者的传播时延
        flood_range_real=flood_B1(1,2)/3e8;
        flood_h(1,1)=t_send_i_j-flood_l(1,1)*(t_send_i_j-flood_delta_t);

        %填充邻居信息矩阵
        flood_neb_history_KF(1,1)=round_KF;
        flood_neb_history_KF(1,2)=t_send_i_j;
        flood_neb_history_KF(1,3)=t_receive_j;
        flood_neb_history_KF(1,4)=t_receive_i;
        neb_list_flood_KF=flood_neb_history_KF;


        %邻居节点第二次接收到该节点消息，进行频率同步与相位同步，更新卡尔曼初始值
    elseif flood_neb_history_KF(1,1)~=0&&flood_neb_history_KF(1,5)==0
        rou_ij=t_receive_j-t_send_i_j;
        rou_ji=t_receive_i-t_send_i_j;
        flood_delta_t=(rou_ji-rou_ij)/2;
        flood_range=(rou_ji+rou_ij)/2;%本次的传播时间

        flood_range_old=(flood_neb_history_KF(1,3)-flood_neb_history_KF(1,2)+...
            flood_neb_history_KF(1,4)-flood_neb_history_KF(1,2))/2;%上一次的传播时间

        flood_alpha_ij=(t_send_i_j-flood_neb_history_KF(1,2))/(t_receive_j-flood_range-...
            (flood_neb_history_KF(1,3)-flood_range_old));%alphai/alphaj,i->j

        %更新频率和相位参数
        flood_l(1,1)=(flood_alpha_ij);
        flood_h(1,1)=flood_l(1,1)*t_send_i_j+flood_h(1,1)-flood_l(1,1)*(t_send_i_j-flood_delta_t);

        %更新邻居的信息矩阵，包括卡尔曼参数
        flood_neb_history_KF(1,1)=round_KF;
        flood_neb_history_KF(1,2)=t_send_i_j;
        flood_neb_history_KF(1,3)=t_receive_j;
        flood_neb_history_KF(1,4)=t_receive_i;
        flood_neb_history_KF(1,5)=(1/flood_alpha_ij-1);%αj/αi-1
        flood_neb_history_KF(1,6)=-flood_delta_t;
        flood_neb_history_KF(1,7)=5e-17;%协方差矩阵参数offset
        flood_neb_history_KF(1,8)=0;
        flood_neb_history_KF(1,9)=0;
        flood_neb_history_KF(1,10)=6.25e-18;%skew
        neb_list_flood_KF=flood_neb_history_KF;


        %邻居超过三次接收到该节点消息，进行卡尔曼频率同步+相位同步计算
    else %flood_neb_history(j,1)~=0&&flood_neb_history(j,5)~=0
        rou_ij=t_receive_j-t_send_i_j;
        rou_ji=t_receive_i-t_send_i_j;
        flood_delta_t=(rou_ji-rou_ij)/2;
        flood_range=(rou_ji+rou_ij)/2;

        %构造上一时刻的协方差矩阵
        flood_P=zeros(2,2);
        flood_P(1,1)=flood_neb_history_KF(1,7);
        flood_P(1,2)=flood_neb_history_KF(1,8);
        flood_P(2,1)=flood_neb_history_KF(1,9);
        flood_P(2,2)=flood_neb_history_KF(1,10);

        %构造上一时刻的状态矩阵
        flood_x=zeros(2,1);
        flood_x(1,1)=flood_neb_history_KF(1,6);
        flood_x(2,1)=flood_neb_history_KF(1,5);%αj/αi-1

        %卡尔曼滤波过程
        [flood_xnew,flood_Pnew]=test_clock_sync_kf(flood_x,...
            flood_P, ...
            t_send_i_j-flood_neb_history_KF(1,2), ...
            t_send_i_j, ...
            t_receive_j-flood_range);


        % flood_range_old=(flood_neb_history(j,3)-flood_neb_history(j,2)+...
        %     flood_neb_history(j,4)-flood_neb_history(j,2))/2;%上一次的传播时间
        %
        % flood_alpha_ij=(t_send_i_j-flood_neb_history(j,2))/(t_receive_j-flood_range-...
        %     (flood_neb_history(j,3)-flood_range_old));%alphai/alphaj
        % %更新频率和相位参数
        % flood_l(n,1)=flood_l(j,1)*(flood_alpha_ij);
        % flood_h(n,1)=flood_l(j,1)*t_send_i_j+flood_h(j,1)-flood_l(n,1)*(t_send_i_j-flood_delta_t);

        %更新频率和相位参数
        flood_l(1,1)=(1/(flood_xnew(2,1)+1));
        flood_h(1,1)=t_send_i_j-flood_l(1,1)*(t_send_i_j+flood_xnew(1,1));

        %更新邻居的信息矩阵，包括卡尔曼参数
        flood_neb_history_KF(1,1)=round_KF;
        flood_neb_history_KF(1,2)=t_send_i_j;
        flood_neb_history_KF(1,3)=t_receive_j;
        flood_neb_history_KF(1,4)=t_receive_i;
        flood_neb_history_KF(1,5)=flood_xnew(2,1);%αj/αi-1
        flood_neb_history_KF(1,6)=flood_xnew(1,1);
        flood_neb_history_KF(1,7)=flood_Pnew(1,1);%协方差矩阵参数
        flood_neb_history_KF(1,8)=flood_Pnew(1,2);
        flood_neb_history_KF(1,9)=flood_Pnew(2,1);
        flood_neb_history_KF(1,10)=flood_Pnew(2,2);
        neb_list_flood_KF=flood_neb_history_KF;
    end
    t_master=t_master+1;
    l_KF(round_KF,1)=flood_l(1,1);
end


skew_flood_KF=l_KF*flood_alpha(2,1)-flood_alpha(1,1);

% --- 3. 绘图 (Y轴对数) ---
semilogy(1:200, abs(skew_flood_KF), '-', 'LineWidth', 1.5);
hold on;
title('卡尔曼滤波跟踪性能');
xlabel('同步轮次');
ylabel('主从频率差绝对值');
grid on;       % 对数坐标下 grid 很有用


%% 最小二乘线性回归滤波

%初始化2节点收到的信息
%每个矩阵元素代表：泛洪轮次，节点上次发送时间，上次自己接收时间，上次对方接收时间,传播时延
neb_list_flood_LR=zeros(200,5);
flood_l=ones(1,1);%泛洪频率补偿
flood_h=zeros(1,1);%泛洪相位补偿
t_master=flood_alpha(1,1)+flood_beta(1,1);%主节点初次泛洪时的本地物理时间
flood_alpha=fre_initial;%泛洪同步频率初始值
flood_beta=pha_initial;%泛洪同步相位初始值

%主节点以1s为周期与邻居进行同步
for round_LR=1:200

   

    %发送消息时刻的距离矩阵
    flood_B1=history.dist_mat{ceil((t_master-flood_beta(1,1))/(flood_alpha(1,1)*time_step))};

    t_send_i_j=t_master+5e-9*randn;%二者约定的发送时间
    %邻居接收时间计算
    t_receive_j=(((t_master-...
        flood_beta(1,1))/(flood_alpha(1,1)))+flood_B1(1,2)/3e8)*flood_alpha(2,1)+...
        flood_beta(2,1)+5e-9*randn;
    %自身接收时间计算
    t_receive_i=(((t_master-...
        flood_beta(2,1))/(flood_alpha(2,1)))+flood_B1(1,2)/3e8)*flood_alpha(1,1)+...
        flood_beta(1,1)+5e-9*randn;

    % 每个矩阵元素代表：泛洪轮次，节点上次发送时间，上次自己接收时间，上次对方接收时间
    flood_neb_history_LR=neb_list_flood_LR;%邻居节点的历史信息矩阵

    %邻居节点第一次接收到该节点信息，只进行相位同步
    if flood_neb_history_LR(1,1)==0
        rou_ij=t_receive_j-t_send_i_j;
        rou_ji=t_receive_i-t_send_i_j;
        flood_delta_t=(rou_ji-rou_ij)/2;%j比i慢多少
        flood_range=(rou_ji+rou_ij)/2;%二者的传播时延
        flood_range_real=flood_B1(1,2)/3e8;
        flood_h(1,1)=t_send_i_j-flood_l(1,1)*(t_send_i_j-flood_delta_t);

        %填充邻居信息矩阵
        flood_neb_history_LR(round_LR,1)=round_LR;
        flood_neb_history_LR(round_LR,2)=t_send_i_j;
        flood_neb_history_LR(round_LR,3)=t_receive_j;
        flood_neb_history_LR(round_LR,4)=t_receive_i;
        flood_neb_history_LR(round_LR,5)=flood_range;
        neb_list_flood_LR=flood_neb_history_LR;


        %邻居节点第二次及以上接收到该节点消息，进行频率同步与相位同步，利用线性回归
    else 
        rou_ij=t_receive_j-t_send_i_j;
        rou_ji=t_receive_i-t_send_i_j;
        flood_delta_t=(rou_ji-rou_ij)/2;
        flood_range=(rou_ji+rou_ij)/2;%本次的传播时间

        i_data=zeros(1,round_LR);
        j_data=zeros(1,round_LR);

        for jj=1:round_LR-1
            i_data(1,jj)=flood_neb_history_LR(jj,2);
            j_data(1,jj)=flood_neb_history_LR(jj,3)-flood_neb_history_LR(jj,5);
        end
        i_data(1,round_LR)=t_send_i_j;
        j_data(1,round_LR)=t_receive_j-flood_range;
        p = polyfit(j_data, i_data, 1);
        flood_alpha_ij=p(1);%alphai/alphaj,i->j

        %更新频率和相位参数
        flood_l(1,1)=(flood_alpha_ij);
        flood_h(1,1)=flood_l(1,1)*t_send_i_j+flood_h(1,1)-flood_l(1,1)*(t_send_i_j-flood_delta_t);

        %更新邻居的信息矩阵，包括卡尔曼参数
        flood_neb_history_LR(round_LR,1)=round_LR;
        flood_neb_history_LR(round_LR,2)=t_send_i_j;
        flood_neb_history_LR(round_LR,3)=t_receive_j;
        flood_neb_history_LR(round_LR,4)=t_receive_i;
        flood_neb_history_LR(round_LR,5)=flood_range;%αj/αi-1
        
        neb_list_flood_LR=flood_neb_history_LR;
    end

       
    t_master=t_master+1;
    l_LR(round_LR,1)=flood_l(1,1);
end


skew_flood_LR=l_LR*flood_alpha(2,1)-flood_alpha(1,1);
%skew_flood_LR=l_LR-alpha_master/flood_alpha(2,1);

l_LR_ij=flood_alpha(1,1)/flood_alpha(2,1);

% --- 3. 绘图 (Y轴对数) ---
semilogy(1:200, abs(skew_flood_LR), '-', 'LineWidth', 1.5);
hold on;
grid on;       % 对数坐标下 grid 很有用


%% 低通滤波

%初始化2节点收到的信息
%每个矩阵元素代表：泛洪轮次，节点上次发送时间，上次自己接收时间，上次对方接收时间，上次的alphaij
neb_list_flood_LP=zeros(1,5);
flood_l=ones(1,1);%泛洪频率补偿
flood_h=zeros(1,1);%泛洪相位补偿
t_master=flood_alpha(1,1)+flood_beta(1,1);%主节点初次泛洪时的本地物理时间
flood_alpha=fre_initial;%泛洪同步频率初始值
flood_beta=pha_initial;%泛洪同步相位初始值

%主节点以1s为周期与邻居进行同步
for round_LP=1:200

   
    

    %发送消息时刻的距离矩阵
    flood_B1=history.dist_mat{ceil((t_master-flood_beta(1,1))/(flood_alpha(1,1)*time_step))};
    t_send_i_j=t_master+5e-9*randn;%二者约定的发送时间
    %邻居接收时间计算
    t_receive_j=(((t_master-...
        flood_beta(1,1))/(flood_alpha(1,1)))+flood_B1(1,2)/3e8)*flood_alpha(2,1)+...
        flood_beta(2,1)+5e-9*randn;
    %自身接收时间计算
    t_receive_i=(((t_master-...
        flood_beta(2,1))/(flood_alpha(2,1)))+flood_B1(1,2)/3e8)*flood_alpha(1,1)+...
        flood_beta(1,1)+5e-9*randn;

    % 每个矩阵元素代表：泛洪轮次，节点上次发送时间，上次自己接收时间，上次对方接收时间，
    % 卡尔曼滤波alpha，卡尔曼滤波beta
    flood_neb_history_LP=neb_list_flood_LP;%邻居节点的历史信息矩阵

    %邻居节点第一次接收到该节点信息，只进行相位同步
    if flood_neb_history_LP(1,1)==0
        rou_ij=t_receive_j-t_send_i_j;
        rou_ji=t_receive_i-t_send_i_j;
        flood_delta_t=(rou_ji-rou_ij)/2;%j比i慢多少
        flood_range=(rou_ji+rou_ij)/2;%二者的传播时延
        flood_range_real=flood_B1(1,2)/3e8;
        flood_h(1,1)=t_send_i_j-flood_l(1,1)*(t_send_i_j-flood_delta_t);

        %填充邻居信息矩阵
        flood_neb_history_LP(1,1)=round_LP;
        flood_neb_history_LP(1,2)=t_send_i_j;
        flood_neb_history_LP(1,3)=t_receive_j;
        flood_neb_history_LP(1,4)=t_receive_i;
        neb_list_flood_LP=flood_neb_history_LP;


        %邻居节点第二次接收到该节点消息，进行频率同步与相位同步，更新alphaij
    elseif flood_neb_history_LP(1,1)~=0&&flood_neb_history_LP(1,5)==0
        rou_ij=t_receive_j-t_send_i_j;
        rou_ji=t_receive_i-t_send_i_j;
        flood_delta_t=(rou_ji-rou_ij)/2;
        flood_range=(rou_ji+rou_ij)/2;%本次的传播时间

        flood_range_old=(flood_neb_history_LP(1,3)-flood_neb_history_LP(1,2)+...
            flood_neb_history_LP(1,4)-flood_neb_history_LP(1,2))/2;%上一次的传播时间

        flood_alpha_ij=(t_send_i_j-flood_neb_history_LP(1,2))/(t_receive_j-flood_range-...
            (flood_neb_history_LP(1,3)-flood_range_old));%alphai/alphaj,i->j

        %更新频率和相位参数
        flood_l(1,1)=(flood_alpha_ij);
        flood_h(1,1)=flood_l(1,1)*t_send_i_j+flood_h(1,1)-flood_l(1,1)*(t_send_i_j-flood_delta_t);

        %更新邻居的信息矩阵，包括卡尔曼参数
        flood_neb_history_LP(1,1)=round_LP;
        flood_neb_history_LP(1,2)=t_send_i_j;
        flood_neb_history_LP(1,3)=t_receive_j;
        flood_neb_history_LP(1,4)=t_receive_i;
        flood_neb_history_LP(1,5)=flood_alpha_ij;
    
        neb_list_flood_LP=flood_neb_history_LP;


        %邻居超过三次接收到该节点消息，进行卡尔曼频率同步+相位同步计算
    else %flood_neb_history(j,1)~=0&&flood_neb_history(j,5)~=0
        rou_ij=t_receive_j-t_send_i_j;
        rou_ji=t_receive_i-t_send_i_j;
        flood_delta_t=(rou_ji-rou_ij)/2;
        flood_range=(rou_ji+rou_ij)/2;

       

        flood_range_old=(flood_neb_history_LP(1,3)-flood_neb_history_LP(1,2)+...
            flood_neb_history_LP(1,4)-flood_neb_history_LP(1,2))/2;%上一次的传播时间

        flood_alpha_ij=(t_send_i_j-flood_neb_history_LP(1,2))/(t_receive_j-flood_range-...
            (flood_neb_history_LP(1,3)-flood_range_old));%alphai/alphaj
        flood_alpha_ij=flood_neb_history_LP(1,5)*0.5+flood_alpha_ij*0.5;
        %更新频率和相位参数
        flood_l(1,1)=(flood_alpha_ij);
        flood_h(1,1)=t_send_i_j-flood_l(1,1)*(t_send_i_j-flood_delta_t);

       

        %更新邻居的信息矩阵，包括卡尔曼参数
        flood_neb_history_LP(1,1)=round_LP;
        flood_neb_history_LP(1,2)=t_send_i_j;
        flood_neb_history_LP(1,3)=t_receive_j;
        flood_neb_history_LP(1,4)=t_receive_i;
        flood_neb_history_LP(1,5)=flood_alpha_ij;%αj/αi-1
        
        neb_list_flood_LP=flood_neb_history_LP;
    end
    t_master=t_master+1;
    l_LP(round_LP,1)=flood_l(1,1);
end

skew_flood_LP=l_LP*flood_alpha(2,1)-flood_alpha(1,1);

% --- 3. 绘图 (Y轴对数) ---
semilogy(1:200, abs(skew_flood_LP), '-', 'LineWidth', 1.5);
hold on;
grid on;       % 对数坐标下 grid 很有用
