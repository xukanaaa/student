
%泛洪伪距时间同步
%分为MTS粗同步，主节点选举，泛洪细同步三个部分
%细同步部分，相位同步用伪距之差，频率同步在多次接收到同一个节点消息基础上，引入卡尔曼滤波，进行降噪与跟踪


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







%=====================================================
flood_alpha=fre_initial;%泛洪同步频率初始值
flood_beta=pha_initial;%泛洪同步相位初始值
master_node=max_node;%主节点编号
begin_round=round_cluster;%开始泛洪的轮次
flood_round=1;%当前泛洪轮次
t_master=(begin_round+1)*flood_alpha(max_node,1)+flood_beta(max_node,1);%主节点初次泛洪时的本地物理时间
flood_l=ones(num_drones,1);%泛洪频率补偿
flood_h=zeros(num_drones,1);%泛洪相位补偿
flood_h_k=zeros(num_drones,1);


%初始化每个节点收到的信息
%每个矩阵元素代表：泛洪轮次，节点上次发送时间，上次自己接收时间，上次对方接收时间，卡尔曼滤波alpha，卡尔曼滤波beta
%最后四列表示协方差矩阵11,12,21,22位置的值
neb_list_flood=cell(num_drones,1);
for i=1:num_drones
    neb_list_flood{i}=zeros(num_drones,10);
end

begin_time=zeros(num_drones,1);%表示每个节点当前泛洪的初始发送时间
begin_time(master_node,1)=t_master;

%初始化节点层级矩阵与轮次矩阵
%初始为50（主节点为0，1）
%一轮泛洪中，从0-50按层级泛洪
%如果邻居轮次小于自身轮次，则邻居层级=自身+1，邻居轮次=自身轮次
%如果邻居轮次等于自身轮次，层级大于自身层级+1，则接收自己消息，否则忽略
%如果邻居轮次大于自身轮次，则忽略
level=zeros(num_drones,2);
for i=1:num_drones
    level(i,1)=50;
    if i==master_node
        level(i,1)=0;
        level(i,2)=1;
    end
end

x_flood=zeros(num_drones,200);
y_flood=zeros(num_drones,200);
delta_max_flood=zeros(200,1);
%主节点以4s为周期与邻居进行同步，每个节点在收到同步消息后间隔0.5s与自己的邻居进行同步
%对一个节点而言，相邻邻居同步间隔为10毫秒
while t_master<=190

    %每一轮从0层次开始逐层往外泛洪
    for i=0:50
        for j=1:num_drones
            %属于当前泛洪轮次且属于当前泛洪层次的节点开始与邻居通信
            if level(j,2)==flood_round&&level(j,1)==i
                t_send=begin_time(j,1);
                %发送消息时刻的连接矩阵
                flood_A=history.adj_mat{ceil((t_send-flood_beta(j,1))/(flood_alpha(j,1)*time_step))};
                %发送消息时刻的距离矩阵
                flood_B1=history.dist_mat{ceil((t_send-flood_beta(j,1))/(flood_alpha(j,1)*time_step))};
                flood_nebcount=0;%已经完成通信的邻居数量，用来确定每个邻居的发送时间
                for n=1:num_drones
                    %与当前发送消息节点有通信链路的邻居n
                    if flood_A(j,n)==1
                        %level(,1):层级，level(,2):轮次
                        %邻居泛洪轮次小于自身泛洪轮次，邻居进行同步
                        %邻居泛洪轮次等于自身泛洪轮次且邻居层级大于自身层级，邻居进行同步
                        if (level(n,2)<level(j,2))||(level(n,2)==level(j,2)&&level(n,1)>level(j,1))
                            t_send_i_j=t_send+0.01*flood_nebcount+5e-9*randn;%二者约定的发送时间
                            %邻居接收时间计算
                            t_receive_j=(((t_send+0.01*flood_nebcount-...
                                flood_beta(j,1))/(flood_alpha(j,1)))+flood_B1(j,n)/3e8)*flood_alpha(n,1)+...
                                flood_beta(n,1)+5e-9*randn;
                            %自身接收时间计算
                            t_receive_i=(((t_send+0.01*flood_nebcount-...
                                flood_beta(n,1))/(flood_alpha(n,1)))+flood_B1(j,n)/3e8)*flood_alpha(j,1)+...
                                flood_beta(j,1)+5e-9*randn;

                            % 每个矩阵元素代表：泛洪轮次，节点上次发送时间，上次自己接收时间，上次对方接收时间，
                            % 卡尔曼滤波alpha，卡尔曼滤波beta
                            flood_neb_history=neb_list_flood{n};%邻居节点的历史信息矩阵

                            %邻居节点第一次接收到该节点信息，只进行相位同步
                            if flood_neb_history(j,1)==0
                                rou_ij=t_receive_j-t_send_i_j;
                                rou_ji=t_receive_i-t_send_i_j;
                                flood_delta_t=(rou_ji-rou_ij)/2;%j比i慢多少
                                flood_range=(rou_ji+rou_ij)/2;%二者的传播时延
                                flood_range_real=flood_B1(n,j)/3e8;
                                flood_h(n,1)=flood_l(j,1)*t_send_i_j+flood_h(j,1)-flood_l(n,1)*(t_send_i_j-flood_delta_t);
                                %flood_h(n,1)=flood_delta_t+flood_h(j,1);
                                flood_h_k(n,1)=flood_l(j,1)*t_send_i_j+flood_h(j,1)-flood_l(n,1)*(t_send_i_j-flood_delta_t);

                                flood_range_ij(j,n)=flood_range;
                                %填充邻居信息矩阵
                                flood_neb_history(j,1)=level(j,2);
                                flood_neb_history(j,2)=t_send_i_j;
                                flood_neb_history(j,3)=t_receive_j;
                                flood_neb_history(j,4)=t_receive_i;
                                neb_list_flood{n}=flood_neb_history;

                                %邻居层级与轮次更新
                                level(n,1)=level(j,1)+1;
                                level(n,2)=level(j,2);


                                %邻居节点第二次接收到该节点消息，进行频率同步与相位同步，更新卡尔曼初始值
                            elseif flood_neb_history(j,1)~=0&&flood_neb_history(j,5)==0
                                rou_ij=t_receive_j-t_send_i_j;
                                rou_ji=t_receive_i-t_send_i_j;
                                flood_delta_t=(rou_ji-rou_ij)/2;
                                flood_range=(rou_ji+rou_ij)/2;%本次的传播时间
                                flood_range_real=flood_B1(n,j)/3e8;

                                flood_range_old=(flood_neb_history(j,3)-flood_neb_history(j,2)+...
                                    flood_neb_history(j,4)-flood_neb_history(j,2))/2;%上一次的传播时间

                                flood_alpha_ij=(t_send_i_j-flood_neb_history(j,2))/(t_receive_j-flood_range-...
                                    (flood_neb_history(j,3)-flood_range_old));%alphai/alphaj

                                %更新频率和相位参数
                                flood_l(n,1)=flood_l(j,1)*(flood_alpha_ij);
                                flood_h(n,1)=flood_l(j,1)*t_send_i_j+flood_h(j,1)-flood_l(n,1)*(t_send_i_j-flood_delta_t);
                                flood_h_k(n,1)=flood_l(j,1)*t_send_i_j+flood_h(j,1)-flood_l(n,1)*(t_send_i_j-flood_delta_t);

                                %更新邻居的信息矩阵，包括卡尔曼参数
                                flood_neb_history(j,1)=level(j,2);
                                flood_neb_history(j,2)=t_send_i_j;
                                flood_neb_history(j,3)=t_receive_j;
                                flood_neb_history(j,4)=t_receive_i;
                                flood_neb_history(j,5)=(1/flood_alpha_ij-1);%αj/αi-1
                                flood_neb_history(j,6)=-flood_delta_t;
                                flood_neb_history(j,7)=5e-17;%协方差矩阵参数offset
                                flood_neb_history(j,8)=0;
                                flood_neb_history(j,9)=0;
                                flood_neb_history(j,10)=6.25e-18;%skew
                                neb_list_flood{n}=flood_neb_history;

                                %更新邻居层级与轮次
                                level(n,1)=level(j,1)+1;
                                level(n,2)=level(j,2);


                                %邻居超过三次接收到该节点消息，进行卡尔曼频率同步+相位同步计算
                            else %flood_neb_history(j,1)~=0&&flood_neb_history(j,5)~=0
                                rou_ij=t_receive_j-t_send_i_j;
                                rou_ji=t_receive_i-t_send_i_j;
                                flood_delta_t=(rou_ji-rou_ij)/2;
                                flood_range=(rou_ji+rou_ij)/2;
                                flood_range_real=flood_B1(n,j)/3e8;

                                %构造上一时刻的协方差矩阵
                                flood_P=zeros(2,2);
                                flood_P(1,1)=flood_neb_history(j,7);
                                flood_P(1,2)=flood_neb_history(j,8);
                                flood_P(2,1)=flood_neb_history(j,9);
                                flood_P(2,2)=flood_neb_history(j,10);

                                %构造上一时刻的状态矩阵
                                flood_x=zeros(2,1);
                                flood_x(1,1)=flood_neb_history(j,6);
                                flood_x(2,1)=flood_neb_history(j,5);%αj/αi-1

                                %卡尔曼滤波过程
                                [flood_xnew,flood_Pnew]=clock_sync_kf(flood_x,...
                                    flood_P, ...
                                    t_send_i_j-flood_neb_history(j,2), ...
                                    t_send_i_j, ...
                                    t_send_i_j-flood_delta_t);


                                % flood_range_old=(flood_neb_history(j,3)-flood_neb_history(j,2)+...
                                %     flood_neb_history(j,4)-flood_neb_history(j,2))/2;%上一次的传播时间
                                % 
                                % flood_alpha_ij=(t_send_i_j-flood_neb_history(j,2))/(t_receive_j-flood_range-...
                                %     (flood_neb_history(j,3)-flood_range_old));%alphai/alphaj
                                % %更新频率和相位参数
                                % flood_l(n,1)=flood_l(j,1)*(flood_alpha_ij);
                                % flood_h(n,1)=flood_l(j,1)*t_send_i_j+flood_h(j,1)-flood_l(n,1)*(t_send_i_j-flood_delta_t);

                                %更新频率和相位参数
                                flood_l(n,1)=(1/(flood_xnew(2,1)+1))*flood_l(j,1);
                                flood_h(n,1)=flood_l(j,1)*t_send_i_j+flood_h(j,1)-flood_l(n,1)*(t_send_i_j+flood_xnew(1,1));
                                flood_h_k(n,1)=flood_l(j,1)*t_send_i_j+flood_h(j,1)-flood_l(n,1)*(t_send_i_j+flood_xnew(1,1));

                                %更新邻居的信息矩阵，包括卡尔曼参数
                                flood_neb_history(j,1)=level(j,2);
                                flood_neb_history(j,2)=t_send_i_j;
                                flood_neb_history(j,3)=t_receive_j;
                                flood_neb_history(j,4)=t_receive_i;
                                flood_neb_history(j,5)=flood_xnew(2,1);%αj/αi-1
                                flood_neb_history(j,6)=flood_xnew(1,1);
                                flood_neb_history(j,7)=flood_Pnew(1,1);%协方差矩阵参数
                                flood_neb_history(j,8)=flood_Pnew(1,2);
                                flood_neb_history(j,9)=flood_Pnew(2,1);
                                flood_neb_history(j,10)=flood_Pnew(2,2);
                                neb_list_flood{n}=flood_neb_history;

                                %更新邻居层级与轮次
                                level(n,1)=level(j,1)+1;
                                level(n,2)=level(j,2);
                            end

                            begin_time(n,1)=t_send+0.01*flood_nebcount+0.5;%更新其邻居下次自己作为主节点发送消息的时间
                            flood_nebcount=flood_nebcount+1;
                        end
                    end
                end
            end
        end
    end
    flood_round=flood_round+1;
    t_master=t_master+4;
    begin_time(master_node,1)=t_master;

    level(master_node,1)=0;
    level(master_node,2)=flood_round;

    %每一轮结束后进行同步性能评估
    %最终调整的逻辑时钟频偏1000100010001
    x_flood(:,4*(flood_round-1)-3)=flood_l.*flood_alpha;
    %最终调整的逻辑时钟相偏
    y_flood(:,4*(flood_round-1)-3)=flood_l.*flood_beta+flood_h;
    %每轮结束后的时钟偏差
    delta_max_flood(4*(flood_round-1)-3,1)=max(x_flood(:,4*(flood_round-1)-3)*...
        t_master+y_flood(:,4*(flood_round-1)-3))-min(x_flood(:,4*(flood_round-1)-3)*...
        t_master+y_flood(:,4*(flood_round-1)-3));



    x1_flood(4*(flood_round-1)-3,1)=max(x_flood(:,4*(flood_round-1)-3))-min(x_flood(:,4*(flood_round-1)-3));


    y1_flood(4*(flood_round-1)-3,1)=max(y_flood(:,4*(flood_round-1)-3))-min(y_flood(:,4*(flood_round-1)-3));


end
fre_temp=flood_l.*flood_alpha;
pha_temp=flood_l.*flood_beta+flood_h;
pha_temp_k=flood_l.*flood_beta+flood_h_k;
aa_fre_cha=max(fre_temp)-min(fre_temp);
aa_pha_cha=max(pha_temp)-min(pha_temp);
aa_pha_k_cha=max(pha_temp_k)-min(pha_temp_k);

skew_flood=x1_flood';
offset_flood=y1_flood';
clock_flood=delta_max_flood';

% --- 2. 提取与拼接skew (同上一步) ---
skew_part1_data = skew_mts(1:round_cluster);
skew_part1_idx = 1:round_cluster;

skew_part2_data = skew_flood(skew_flood ~= 0); 
count_non_zero = length(skew_part2_data);
skew_part2_idx = 1 : 4 : (1 + (count_non_zero - 1) * 4);

skew_Y_combined = [skew_part1_data, skew_part2_data];
skew_X_combined = [skew_part1_idx, skew_part2_idx+round_cluster]; 

% --- 2. 提取与拼接offset (同上一步) ---
offset_part1_data = offset_mts(1:round_cluster);
offset_part1_idx = 1:round_cluster;

offset_part2_data = offset_flood(offset_flood ~= 0); 
count_non_zero = length(offset_part2_data);
offset_part2_idx = 1 : 4 : (1 + (count_non_zero - 1) * 4);

offset_Y_combined = [offset_part1_data, offset_part2_data];
offset_X_combined = [offset_part1_idx, offset_part2_idx+round_cluster];

% --- 2. 提取与拼接clock (同上一步) ---
clock_part1_data = clock_mts(1:round_cluster);
clock_part1_idx = 1:round_cluster;

clock_part2_data = clock_flood(clock_flood ~= 0); 
count_non_zero = length(clock_part2_data);
clock_part2_idx = 1 : 4 : (1 + (count_non_zero - 1) * 4);

clock_Y_combined = [clock_part1_data, clock_part2_data];
clock_X_combined = [clock_part1_idx, clock_part2_idx+round_cluster];

% --- 3. 绘图 (Y轴对数) ---
semilogy(skew_X_combined, skew_Y_combined, '-o', 'LineWidth', 1.5);
hold on;
% semilogy(offset_X_combined, offset_Y_combined, '-o', 'LineWidth', 1.5);
% hold on;
semilogy(clock_X_combined, clock_Y_combined, '-o', 'LineWidth', 1.5);
hold on;

title('拼接数据 - 对数坐标 (Y-Log Scale)');
xlabel('原始矩阵索引 (Original Index)');
ylabel('数值 (Log Scale)');
grid on;       % 对数坐标下 grid 很有用