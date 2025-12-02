
%% 超广域战术蜂群仿真 (Tactical Scale Swarm)
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
%1-5:5，6-15：10,16-30:15,31-50：20

t_period=1;%发送消息伪周期
time_step = 0.1; % 时间步长（s）
simulation_k=180;%仿真轮数
%不同轮次调整后的频偏和相偏
x=zeros(num_drones,simulation_k);
y=zeros(num_drones,simulation_k);
z=zeros(num_drones,simulation_k);
deltamax=zeros(simulation_k,1);%每轮迭代后全网最大钟差
var_clock=zeros(simulation_k,1);%每轮迭代后全网时钟值的方差
t=zeros(simulation_k,1);%每轮迭代后用于测评的真实时间，此处选为最后一个发送消息的节点每次发送完的时间
round=3;%蒙特卡洛仿真轮数
skew=zeros(round,simulation_k);%未平均之前的最大频偏差值矩阵
skewFinal=zeros(1,simulation_k);%平均后的每一轮频偏差值。
phySkew=80e-6;%物理时钟频偏最大值
phyOffset=1;%物理时钟相位偏移最大值
offset=zeros(round,simulation_k);
clock_var=zeros(round,simulation_k);
offsetFinal=zeros(1,simulation_k);
clock=zeros(round,simulation_k);
clockFinal=zeros(1,simulation_k);


%初始化物理时钟频偏
alpha=1-phySkew+(2*phySkew)*rand(num_drones,1);
%初始化物理时钟相偏
beta=-phyOffset+2*phyOffset*rand(num_drones,1);
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
Neb_list_fre=cell(num_drones,1);
for i=1:num_drones
    Neb_list_fre{i}=zeros(num_drones*2,5);
end

%初始化每个节点本地的相位邻居信息列表(消息轮，(当前邻居物理，当前自身物理)，邻居逻辑时钟频率调整l，邻居逻辑时钟相位调整h
%,)
%前n行表示首次接收到某节点的消息数据
Neb_list_pha=cell(num_drones,1);
for i=1:num_drones
    Neb_list_pha{i}=zeros(num_drones,5);
end

%初始化每个节点收到的头结点得分汇总
Neb_list_score=cell(num_drones,1);
for i=1:num_drones
    Neb_list_score{i}=zeros(num_drones,1);
end

round_cluster=0;

for j=1:simulation_k
    for i=1:num_drones
        [~,d]=sort(t_global_total(:,j));

        %第一轮不更新，从第二轮开始
        if j>=2
            M=MTS_l_h(Neb_list_fre{d(i)},Neb_list_pha{d(i)},j,num_drones,l(d(i),1),h(d(i),1),d(i));%更新逻辑时钟频偏和相偏参数
            %返回值包括频偏补偿，相偏补偿，当前同步状态（0,0.5,1）
            l(d(i),1)=M(1,1);
            h(d(i),1)=M(1,2);
            state_self=Neb_list_fre{d(i)};
            state_self(d(i),5)=M(1,3);
            Neb_list_fre{d(i)}=state_self;

        end
        A=history.adj_mat{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的连接矩阵
        B1=history.dist_mat{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的距离矩阵

        %更新完后对当前节点局部区域同步状态进行判断，如果自己连续两轮未改变且当前轮次收到的邻居也是如此，状态为1，否则为0
        F1=Neb_list_fre{d(i)};%广播节点的自身的消息矩阵
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
        S_self=Neb_list_score{d(i)};
        if S_self(d(i),1)==0&&syn_state==1
            score_self=get_score(d(i),A,num_drones,j);
            S_self(d(i),1)=score_self;
        end

        %计算簇头选举出来时的轮次
        if all(S_self(:))
            if S_self(d(i),1)==min(S_self)&&round_cluster==0
                round_cluster=j;%轮次
                A_cluster=A;%连接矩阵
                B_cluster=B1;%距离矩阵
                S_self_temp=S_self;
            end
        end

        %检查d(i)节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
        for k=1:num_drones
            if A(d(i),k)==1

                F=Neb_list_fre{k};%频率矩阵
                P=Neb_list_pha{k};%相位矩阵
                S=Neb_list_score{k};%邻居的分数矩阵

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

                    Neb_list_fre{k}=F;
                    Neb_list_pha{k}=P;
                    Neb_list_score{k}=S;
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

                    F_self=Neb_list_fre{d(i)};
                    F(d(i)+num_drones,5)=F_self(d(i),5);

                    %更新合并自己的分数矩阵和邻居的分数矩阵
                    S_merge=merge(S,S_self,num_drones);
                    S=S_merge;

                    Neb_list_fre{k}=F;
                    Neb_list_pha{k}=P;
                    Neb_list_score{k}=S;
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
    deltamax(k,1)=max(x(:,k)*t(k,1)+y(:,k))-min(x(:,k)*t(k,1)+y(:,k));
end

%每轮结束后同步后时钟的方差
for k=1:simulation_k
    var_clock(k,1)=var(x(:,k)*t(k,1)+y(:,k));
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
%一共进行5次蒙特卡洛仿真，r=5，每一行表示各轮次同步后的全网最大频率差
skew(1,:)=x1';
offset(1,:)=y1';
clock(1,:)=deltamax';
clock_var(1,:)=var_clock';

%五次仿真求平均值
skewFinal=skew(1,:);
offsetFinal=offset(1,:);
clockFinal=clock(1,:);
clock_varFinal=clock_var(1,:);


format long;
semilogy(1:5:simulation_k,skewFinal(1:5:end),'-s','LineWidth',1.5);
 hold on;
%plot(1:simulation_k,offsetFinal(1:end),'-*','LineWidth',1);
% hold on;
%semilogy(1:simulation_k,clockFinal(1:end),'-*','LineWidth',1);
 semilogy(1:5:simulation_k,clockFinal(1:5:end),'-*','LineWidth',1.5);
xlabel("同步轮数");
ylabel("全网最大时钟差");

% semilogy(1:simulation_k,clock_varFinal(1:end),'-*','LineWidth',1);
grid on;
hold on;


%% 2. 构建图
% 融合 A 和 B：只保留 A 中有连接的边，并赋予 B 的权重
W = double(A_cluster) .* B_cluster;
G = graph(W, 'upper');

%% 3. 绘图 (基础设置)
figure('Color', 'w', 'Position', [100, 100, 900, 700]);

% 使用力导向布局，让距离近的节点聚在一起
h = plot(G, 'Layout', 'force'); 

title('网络拓扑图 (仅显示 1-5 号节点编号)');
axis off;

% --- 样式微调：让背景节点淡化 ---
h.NodeColor = [0.8 0.8 0.8]; % 灰色节点
h.MarkerSize = 5;            % 小节点
h.EdgeColor = [0.9 0.9 0.9]; % 非常淡的边线
h.EdgeAlpha = 0.6;

%% 4. 关键步骤：只显示 1-5 号的编号
% 默认情况下 MATLAB 会显示所有编号。我们需要自定义标签列表。
% --- 修正开始 ---
% 1. 创建一个全是空字符串 "" 的字符串向量
myLabels = strings(num_drones, 1); 

% 2. 填入前 5 个
for i = 1:5
    myLabels(i) = string(i);
end

% 3. 赋值
h.NodeLabel = myLabels;
% --- 修正结束 ---
% 其余位置保持为空 cell，即不显示标签

% 将自定义标签赋值给图对象
h.NodeLabel = myLabels;

%% 5. 高亮 1-5 号节点 (红色 + 大字体)
target_nodes = 1:5;

highlight(h, target_nodes, ...
    'NodeColor', 'r', ...         % 节点变红
    'MarkerSize', 10, ...         % 节点变大
    'NodeLabelColor', 'k', ...    % 编号文字设为黑色 (对比度高)
    'NodeFontSize', 14, ...       % 编号字体加大
    'NodeFontWeight', 'bold');    % 编号字体加粗

disp('绘图完成。仅 1-5 号节点显示编号。');
