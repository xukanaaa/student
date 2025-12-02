%使用首轮和最新轮次消息用来同步，并且对相对运动产生的位置偏移进行补偿，利用线性回归最小二乘滤波
%在粗同步之后尝试利用强化学习来自适应同步周期。
%粗同步阶段有虚拟链接，强化学习阶段不采用虚拟链接
result_zuihou=zeros(3,1);

num_drones = 100; % 无人机数量
space_size = 3000; % 空间大小（长、宽、高均为1km）
min_speed =80; % 最小速度（m/s）
max_speed =100; % 最大速度（m/s）
comm_range = 1500; % 通信范围（m）
time_step = 0.1; % 时间步长（s）
speedChange=2;%无人机每秒速度的该变量（正态分布）
simulation_time = 500; % 模拟时间（s）
t_period=2;%发送消息伪周期
simulation_k=200;%仿真轮数
broad=5;%多重广播一次的广播次数
dt=20e-3;%多重广播间隔
%不同轮次调整后的频偏和相偏
x=zeros(num_drones,simulation_k);
y=zeros(num_drones,simulation_k);
z=zeros(num_drones,simulation_k);
deltamax=zeros(simulation_k,1);%每轮迭代后全网最大钟差
var_clock=zeros(simulation_k,1);%每轮迭代后全网时钟值的方差
t=zeros(simulation_k,1);%每轮迭代后用于测评的真实时间，此处选为最后一个发送消息的节点每次发送完的时间
speed_change=zeros(num_drones,3);%无人机初始化速度扰动
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
used_num=0;%每轮同步使用的数据量
M=zeros(num_drones+1,2);
Var_change=zeros(num_drones,1);%存储上一轮的方差，用来决定同步趋势
State_old=zeros(num_drones,1);%存储上一轮的状态，将（方差，趋势）映射到唯一状态索引
alpha_RL=0.95;%强化学习的学习率
gamma_RL=0.9;%强化学习未来衰减率
Qtable_trance=cell(400,1);%追踪Qtable迭代过程
Q_count=1;
clock_max_eachround=zeros(simulation_k,1);%用来统计强化学习中和原本类似时间网络的最大时钟差
%complexity_RL=zeros(30,1);%强化学习复杂度，每50s对应时间为一块，采用累加
%complexity_RL(1,1)=2500;
complexity_DPTS=[2500;5000;7500;10000;12500;15000;17500;20000;22500;25000;27500;30000;32500;35000;37500;40000;42500;45000;47500;50000];

for r=1:3
    complexity_RL=zeros(30,1);%强化学习复杂度，每50s对应时间为一块，采用累加
    complexity_RL(1,1)=2500;
    % 初始化无人机位置
    drone_positions = space_size * rand(num_drones, 3);
    % 初始化无人机速度
    drone_speeds=zeros(num_drones,3);
    for i=1:num_drones
        for j=1:3
            positive_value=unifrnd(min_speed,max_speed);
            sign=randi([0,1])*2-1;
            drone_speeds(i,j)=sign*positive_value;
        end
    end

    %初始化物理时钟频偏
    alpha=1-phySkew+(2*phySkew)*rand(num_drones,1);
    %初始化物理时钟相偏
    beta=-phyOffset+2*phyOffset*rand(num_drones,1);
    %初始化逻辑时钟频率调整参数
    l=ones(num_drones,1);
    %初始化逻辑时钟相位调整参数
    h=zeros(num_drones,1);
    % 初始化无人机通信邻接矩阵
    comm_matrix = zeros(num_drones);

    %初始化每个离散时间点的距离矩阵
    D=cell(simulation_time/time_step,1);
    for i=1:(simulation_time/time_step)
        D{i}=zeros(num_drones,num_drones);
    end

    %初始化每个离散时间点的连接矩阵
    I=cell(simulation_time/time_step,1);
    for i=1:(simulation_time/time_step)
        I{i}=zeros(num_drones,num_drones);
    end

    % 模拟无人机运动
    for k=1:(simulation_time/time_step)

        %每隔一秒更新无人机速度
        if mod(k,1/time_step)==0
            speed_change=speedChange*randn(num_drones,3);%无人机速度扰动
            drone_speeds=drone_speeds+speed_change;%更新扰动后的速度
        end

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
    end

    t_local_total=zeros(num_drones,simulation_k);
    t_global_total=zeros(num_drones,simulation_k);
    t_global=(t_period*0.9)*sort(rand(num_drones,1));%初始发送全球时间，升序,预留0.1倍周期作为保护jiange
    t_local=alpha.*t_global+beta;%初始发送节点时间
    t_local_total(:,1)=t_local;
    t_global_total(:,1)=t_global;

    %所有轮次中各节点发送消息的本地物理时间
    for m=2:simulation_k
        t_local_total(:,m)=t_local+t_period*(m-1)*alpha;
    end

    %所有轮次中各节点发送消息的全局时间
    for m=2:simulation_k
        t_global_total(:,m)=(t_local_total(:,m)-beta)./alpha;
    end

    %初始化每个节点本地的邻居频率调整信息列表(消息轮，(当前邻居物理，当前自身物理)*广播轮数，当前邻居逻辑时钟频率调整l,
    %flag位，只占该列第一行，下面都是0
    %最后n列：对应邻居节点的二跳邻居信息，每次更新完要清零)
    %前n行表示首次接收到某节点的消息数据
    %后n行表示最新接收到某节点的消息数据
    Neb_list_fre=cell(num_drones,1);
    for i=1:num_drones
        Neb_list_fre{i}=zeros(num_drones*2,broad*2+2+1+num_drones);
    end

    %初始化每个节点本地的相位邻居信息列表(消息轮，(当前邻居物理，当前自身物理)*广播轮数，邻居逻辑时钟频率调整l，邻居逻辑时钟相位调整h
    %,与这个邻居的距离产生的时延（随着运动会实时补偿变化）,最后i列表示这个邻居所带来的虚拟链接的时钟相位差值)
    %前n行表示首次接收到某节点的消息数据
    Neb_list_pha=cell(num_drones,1);
    for i=1:num_drones
        Neb_list_pha{i}=zeros(num_drones,1+broad*2+2+1+num_drones);
    end

    for j=1:25 %前25轮完成初步同步
        for i=1:num_drones
            [~,d]=sort(t_global_total(:,j));


            A=I{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的连接矩阵
            B1=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的距离矩阵
            %B2=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))+1};

            if j>=2
                Vara=zeros(num_drones,1);
                Atemp=Neb_list_pha{d(i)};
                Btemp=zeros(2,1);%第一行表示邻居的逻辑时钟，第二行表示自身的逻辑时钟
                for n=1:num_drones
                    if Atemp(n,1)>=j-1

                        Atemp(n,2)=Atemp(n,2)-Atemp(n,14);%传播时延补偿


                        Btemp(1,1)=Atemp(n,2)*Atemp(n,12)+Atemp(n,13);
                        Btemp(2,1)=Atemp(n,3)*l(d(i),1)+h(d(i),1);
                        Vara(n,1)=Btemp(1,1)-Btemp(2,1);

                    end
                end
                Var_change(d(i),1)=var(Vara(Vara ~= 0));%方差初始化

                %节点状态初始化
                State_old(d(i),1)=state_to_index(Var_change(d(i),1),3);
            end

            %第一轮不更新，从第二轮开始
            if j>=2
                M=distributed_all_l(Neb_list_fre{d(i)},j,num_drones,l(d(i),1));%更新逻辑时钟频偏参数
                %这个返回的M矩阵其实包含了在这一轮更新中收到真实邻居的αij*l的值，后面广播可以放到该节点携带的虚拟信息列表中。
                l(d(i),1)=M(1,1);%更新逻辑时钟频偏参数
                used_num=M(1,2);%取出更新使用的数据量，决定下一轮是否广播求助信息，让接收到该信息的邻居节点广播二跳信息。
                %更新完之后，把右侧n列关于二跳节点的信息清零(频率)
                E=Neb_list_fre{d(i)};
                E(:, 14:end)=0;
                Neb_list_fre{d(i)}=E;

                N=distributed_all_h(Neb_list_pha{d(i)},j,l(d(i),1),h(d(i),1),num_drones);%更新逻辑时钟相偏参数
                h(d(i),1)=N(1,1);%更新逻辑时钟相偏参数
                %更新完之后，把右侧n列关于二跳节点的信息清零(相位)
                E=Neb_list_pha{d(i)};
                E(:, 15:end)=0;
                Neb_list_pha{d(i)}=E;
            end


            F1=Neb_list_fre{d(i)};%广播节点的自身的消息矩阵

            %检查d(i)节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
            for k=1:num_drones
                if A(d(i),k)==1

                    F=Neb_list_fre{k};%频率矩阵
                    P=Neb_list_pha{k};%相位矩阵
                    %对于自己上一轮用于同步邻居小于3时，向所有邻居发送求助，修改邻居的flag为1，
                    %该邻居检测到flag为1时，则在下次广播中附加自己的二跳信息
                    if j>=2 && used_num<=5
                        F(1,13)=1;
                    end
                    %首次接收到的消息（不管是第几轮），存放在前n行，n表示无人机总数
                    if F(d(i),1)==0
                        F(d(i),1)=j;
                        P(d(i),1)=j;
                        %填充邻居列表中首次接收到该节点信息的数据(频率和相位矩阵都有)
                        for k1=1:5
                            F(d(i),2*k1)=t_local_total(d(i),j)+dt*(k1-1)+5e-9*randn;
                            P(d(i),2*k1)=F(d(i),2*k1);
                        end

                        for k2=1:5
                            F(d(i),2*k2+1)=(((t_local_total(d(i),j)+(k2-1)*dt-...
                                beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;
                            P(d(i),2*k2+1)=F(d(i),2*k2+1);
                        end
                        F(d(i),12)=l(d(i),1);
                        P(d(i),12)=l(d(i),1);
                        P(d(i),13)=h(d(i),1);
                        %首次收到相位消息，利用双向消息交换计算距离时延
                        T1=F(d(i),2);
                        T2=F(d(i),3);
                        T3=T2+20e-3;
                        T4=(((T3-beta(k,1))/(alpha(k,1)))+B1(d(i),k)/3e8)*alpha(d(i),1)+beta(d(i),1)+5e-9*randn;
                        Td=(T2+T4-T3-T1)/2;
                        P(d(i),14)=Td;

                        %判断是否往邻居消息列表中加入当前广播节点的二跳信息(频率和相位)，如果广播节点flag=1，则需要添加
                        if F1(1,13)==1
                            for k3=2:num_drones+1
                                if M(k3,1)>=j-1
                                    F(d(i),13+k3-1)=M(k3,2);
                                end
                                if N(k3,1)>=j-1
                                    P(d(i),13+k3)=N(k3,2);
                                end
                            end
                        end
                        Neb_list_fre{k}=F;
                        Neb_list_pha{k}=P;
                    end

                    %非首次最新接受到的消息（不管是那一轮），
                    %频率消息存放在第n+1到2n行，用于计算的数据是首次数据和最新数据
                    %相位消息首先经过Td补偿，然后更新第一轮
                    if F(d(i),1)~=0
                        F(d(i)+num_drones,1)=j;
                        P(d(i),1)=j;
                        for k1=1:5
                            F(d(i)+num_drones,2*k1)=t_local_total(d(i),j)+dt*(k1-1)+5e-9*randn;
                        end

                        for k2=1:5
                            F(d(i)+num_drones,2*k2+1)=(((t_local_total(d(i),j)+(k2-1)*dt-...
                                beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;
                        end
                        deltaTd=compensatePha(F,P,d(i),num_drones);%对新消息的相位消息的d进行动态更新
                        %更新相位矩阵的值
                        for k1=1:5
                            P(d(i),2*k1)=F(d(i)+num_drones,2*k1);
                            P(d(i),2*k1+1)=F(d(i)+num_drones,2*k1+1);
                        end
                        P(d(i),12)=l(d(i),1);
                        P(d(i),13)=h(d(i),1);
                        P(d(i),14)=P(d(i),14)+deltaTd;
                        F(d(i)+num_drones,12)=l(d(i),1);

                        %判断是否往邻居消息列表中加入当前广播节点的二跳信息(频率和相位)，如果广播节点flag=1，则需要添加
                        if F1(1,13)==1
                            for k3=2:num_drones+1
                                if M(k3,1)>=j-1
                                    F(d(i)+num_drones,13+k3-1)=M(k3,2);
                                end
                                if N(k3,1)>=j-1
                                    P(d(i),13+k3)=N(k3,2);
                                end
                            end
                        end
                        Neb_list_fre{k}=F;
                        Neb_list_pha{k}=P;
                    end
                end
            end
            %当所有的消息都广播完毕后，广播节点的flag要重新置零
            F1(1,13)=0;
            Neb_list_fre{d(i)}=F1;
        end
        %最终调整的逻辑时钟频偏
        x(:,j)=l.*alpha;
        %最终调整的逻辑时钟相偏
        y(:,j)=l.*beta+h;
        z(:,j)=h;
    end

    for i=26:simulation_k
        x(:,i)=l.*alpha;
    end
    for i=26:simulation_k
        y(:,i)=l.*beta+h;
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
    skew(r,:)=x1';
    offset(r,:)=y1';
    clock(r,:)=deltamax';
    clock_var(r,:)=var_clock';


    %初步同步之后，采用强化学习方法自适应调整周期。权衡同步精度和通信开销
    %这种情况下，频率同步和相位同步无需使用虚拟链接
    %不采用同步轮次，每个节点作为独立智能体学习不同情况下的周期调整策略。
    %动作空间为-5，-2，-1,0,1,2,5，动作空间大小为7
    %方差离散点为1~99*10^-11,一共50个状态，其中最后25个状态方差已经不可接受，每个状态对应up，down，equal三个趋势，状态空间大小为150

    Qtable=zeros(150,7);%初始化Q表格
    Tsend=zeros(num_drones,3);%存储每个节点下次广播的全局时间、调整完的同步周期和调整周期的动作
    Tchange=100*ones(num_drones,2000);%存储每个节点周期的变化调整过程。
    Sold=zeros(num_drones,2);%存储上一个状态对（方差，趋势）


    %初始化每个节点本地的邻居频率调整信息列表(消息轮，(当前邻居物理，当前自身物理)*广播轮数，当前邻居逻辑时钟频率调整l
    %前n行表示首次接收到某节点的消息数据
    %后n行表示最新接收到某节点的消息数据
    Neb_list_fre_RL=cell(num_drones,1);
    for i=1:num_drones
        Neb_list_fre_RL{i}=zeros(num_drones*2,broad*2+2);
    end

    %初始化每个节点本地的相位邻居信息列表(消息轮，(当前邻居物理，当前自身物理)*广播轮数，邻居逻辑时钟频率调整l，邻居逻辑时钟相位调整h
    %,与这个邻居的距离产生的时延（随着运动会实时补偿变化）)
    %前n行表示首次接收到某节点的消息数据
    Neb_list_pha_RL=cell(num_drones,1);
    for i=1:num_drones
        Neb_list_pha_RL{i}=zeros(num_drones,1+broad*2+2+1);
    end
    %每次从中选取最小值节点广播（广播前进行同步），广播过程修改直连邻居参数，广播后修改下次广播时间
    %每次进行同步时，通过比较本次与上次同步前方差，获得上次行为后系统状态，利用奖励函数更新上一状态上一动作的Q值

    %初始化发射时间Tsend，将其赋值为粗同步之后的下一轮全局时间，并且当前周期设置为2，下次调整设置为0
    Tsend(:,1)=(t_local_total(:,26)-beta)./alpha;
    Tsend(:,2)=2 * ones(num_drones, 1);
    Tsend(:,3)=zeros(num_drones, 1);

    %开始Qlearning算法迭代，每次选取当前发射时间最小值节点，先广播，再更新
    while  min(Tsend(:, 1))<500

        %先取广播时间最小节点广播，广播完成后按照2,3列数据更新下次广播全局时间
        [min_value, min_index] = min(Tsend(:, 1));
        %获取广播时刻的连接矩阵与距离矩阵
        A_RL=I{ceil(min_value/time_step)};
        B1_RL=D{ceil(min_value/time_step)};

        %更新邻居的频率同步矩阵和相位同步矩阵
        %频率调整第一列为flag位，1表示最新消息，每次同步后要把为1的flag置零。
        %检查min_index节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
        for k=1:num_drones
            if A_RL(min_index,k)==1

                F_RL=Neb_list_fre_RL{k};%频率矩阵
                P_RL=Neb_list_pha_RL{k};%相位矩阵

                %首次接收到的消息（不管是第几轮），存放在前n行，n表示无人机总数
                if F_RL(min_index,2)==0
                    F_RL(min_index,1)=1;
                    P_RL(min_index,1)=1;
                    %填充邻居列表中首次接收到该节点信息的数据(频率和相位矩阵都有)
                    for k1=1:5
                        F_RL(min_index,2*k1)=min_value*alpha(min_index,1)+beta(min_index,1)+dt*(k1-1)+5e-9*randn;
                        P_RL(min_index,2*k1)=F_RL(min_index,2*k1);
                    end

                    for k2=1:5
                        F_RL(min_index,2*k2+1)=(((min_value*alpha(min_index,1)+beta(min_index,1)+(k2-1)*dt-...
                            beta(min_index,1))/(alpha(min_index,1)))+B1_RL(min_index,k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;
                        P_RL(min_index,2*k2+1)=F_RL(min_index,2*k2+1);
                    end
                    F_RL(min_index,12)=l(min_index,1);
                    P_RL(min_index,12)=l(min_index,1);
                    P_RL(min_index,13)=h(min_index,1);

                    %首次收到相位消息，利用双向消息交换计算距离时延
                    T1=F_RL(min_index,2);
                    T2=F_RL(min_index,3);
                    T3=T2+20e-3;
                    T4=(((T3-beta(k,1))/(alpha(k,1)))+B1_RL(min_index,k)/3e8)*alpha(min_index,1)+beta(min_index,1)+5e-9*randn;
                    Td=(T2+T4-T3-T1)/2;
                    P_RL(min_index,14)=Td;

                    Neb_list_fre_RL{k}=F_RL;
                    Neb_list_pha_RL{k}=P_RL;
                end

                %非首次最新接受到的消息（不管是那一轮），
                %频率消息存放在第n+1到2n行，用于计算的数据是首次数据和最新数据
                %相位消息首先经过Td补偿，然后更新第一轮
                if F_RL(min_index,2)~=0
                    F_RL(min_index+num_drones,1)=1;
                    P_RL(min_index,1)=1;
                    for k1=1:5
                        F_RL(min_index+num_drones,2*k1)=min_value*alpha(min_index,1)+beta(min_index,1)+dt*(k1-1)+5e-9*randn;
                    end

                    for k2=1:5
                        F_RL(min_index+num_drones,2*k2+1)=(((min_value*alpha(min_index,1)+beta(min_index,1)+(k2-1)*dt-...
                            beta(min_index,1))/(alpha(min_index,1)))+B1_RL(min_index,k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;
                    end
                    deltaTd_RL=compensatePha_RL(F_RL,P_RL,min_index,num_drones);%对新消息的相位消息的d进行动态更新
                    %更新相位矩阵的值
                    for k1=1:5
                        P_RL(min_index,2*k1)=F_RL(min_index+num_drones,2*k1);
                        P_RL(min_index,2*k1+1)=F_RL(min_index+num_drones,2*k1+1);
                    end
                    P_RL(min_index,12)=l(min_index,1);
                    P_RL(min_index,13)=h(min_index,1);
                    P_RL(min_index,14)=P_RL(min_index,14)+deltaTd_RL;
                    F_RL(min_index+num_drones,12)=l(min_index,1);

                    Neb_list_fre_RL{k}=F_RL;
                    Neb_list_pha_RL{k}=P_RL;
                end
            end
        end

        %更新Tchange矩阵,用来追踪每个节点每次的动作
        for k3=1:200
            if (Tchange(min_index,k3)==100)
                Tchange(min_index,k3)=Tsend(min_index,3);
                break;
            end
        end

        %第一步：计算方差，获取当前状态，获取上一个动作的奖励，更新上一个状态的Qtable
        P_var_RL=Neb_list_pha_RL{min_index};
        if any(P_var_RL(:) ~= 0) && get_var(P_var_RL,num_drones,l(min_index,1),h(min_index,1))~=0
            var_curr_RL=get_var(P_var_RL,num_drones,l(min_index,1),h(min_index,1));%从相位同步矩阵获取前段时间的时钟值方差，作为当前状态
        else
            var_curr_RL=Var_change(min_index,1);
        end
        var_old_RL=Var_change(min_index,1);%上一轮的方差
        state_old_RL=State_old(min_index,1);%上一轮的转态索引
        %state_curr_RL=0;%这一轮的状态索引初始值

        %获取当前状态索引
        if var_curr_RL-var_old_RL>5e-11 %方差变大
            state_curr_RL=state_to_index(var_curr_RL,1);%状态转索引函数
        elseif var_old_RL-var_curr_RL>5e-11 %方差变小
            state_curr_RL=state_to_index(var_curr_RL,2);
        else %方差几乎保持不变
            state_curr_RL=state_to_index(var_curr_RL,3);
        end

        %将上一轮的方差和状态索引更新为此次的，方便下次使用
        Var_change(min_index,1)=var_curr_RL;
        State_old(min_index,1)=state_curr_RL;

        %从精度，周期，趋势三个方面计算上一个状态采取的动作的奖励R
        R_RL=get_r(Tsend(min_index,3),var_curr_RL,state_curr_RL);

        %Qlearning中的TD方式更新前一个状态的QTable值
        action_old_RL=action_to_index(Tsend(min_index,3));%上一个动作转索引
        max_Q=max(Qtable(state_curr_RL, :));
        Qtable(state_old_RL,action_old_RL)=Qtable(state_old_RL,action_old_RL)+...
            alpha_RL*(R_RL+gamma_RL*max_Q-Qtable(state_old_RL,action_old_RL));

        %第二步：根据贪婪策略选取下次动作，并且更新Tsend的动作，周期和下次发射时间
        action_curr_RL=greedy(state_curr_RL,Qtable);
        Tsend(min_index,3)=action_curr_RL;
        %按照Tsend3列的值更新第一列下次的广播时间和周期
        if Tsend(min_index,2)+Tsend(min_index,3)>2 %周期缩短没有越界
            Tsend(min_index,2)=Tsend(min_index,2)+Tsend(min_index,3);
            Tsend(min_index,1)=(Tsend(min_index,1)*alpha(min_index,1)+beta(min_index,1)+...
                Tsend(min_index,2)-beta(min_index,1))/alpha(min_index,1);
        else
            Tsend(min_index,2)=2;
            Tsend(min_index,1)=(Tsend(min_index,1)*alpha(min_index,1)+beta(min_index,1)+...
                Tsend(min_index,2)-beta(min_index,1))/alpha(min_index,1);
        end

        Qtable_trance{Q_count}=Qtable;
        Q_count=Q_count+1;

        %第三步：完成广播之后，利用自身的频率和相位调整矩阵进行同步
        F_self_RL=Neb_list_fre_RL{min_index};
        P_self_RL=Neb_list_pha_RL{min_index};
        if any(F_self_RL(:) ~= 0) %频率矩阵非空才考虑进行同步
            [l_change,F_change]=distributed_all_l_RL(F_self_RL,num_drones,l(min_index,1));%更新逻辑时钟频偏参数
            %第一个参数表示频偏调整值，第二个表示同步之后的邻居矩阵（要把已经用过的信息标记）
            l(min_index,1)=(1-0.05-0.985^(fix(min_value/2)))*l(min_index,1)+(0.05+0.985^(fix(min_value/2)))*l_change;
            Neb_list_fre_RL{min_index}=F_change;

            [h_change,P_change]=distributed_all_h_RL(P_self_RL,l(min_index,1),h(min_index,1),num_drones);%更新逻辑时钟相偏参数
            h(min_index,1)=(1-0.05-0.985^(fix(min_value/2)))*h(min_index,1)+(0.05+0.985^(fix(min_value/2)))*h_change;%更新逻辑时钟相偏参数
            Neb_list_pha_RL{min_index}=P_change;
        end

        %这个变量其实看的是400s左右网络的最大时钟差
        if min(Tsend(:, 1))>390&&min(Tsend(:, 1))<410
            global_temp=min(Tsend(:, 1));
            local_temp=global_temp*alpha+beta;
            compa_temp=local_temp.*l+h;
            result_temp=max(compa_temp) - min(compa_temp);
            result_zuihou(r,1)=result_temp;
        end
        
        %计算强化学习对应时间点的最大时钟差
        round_send=ceil(min_value/2);%当前同步时刻对应的强化学习之前的同步轮
        if round_send<=200
            local_time_RL=min_value*alpha+beta;
            compensate_time_RL=local_time_RL.*l+h;
            delta_clock_RL=max(compensate_time_RL)-min(compensate_time_RL);
            if(clock_max_eachround(round_send,1)==0 )
                clock_max_eachround(round_send,1)=delta_clock_RL;
            end
        end

        %计算强化学习通信开销
        if min_value<=1000
            complexity_index=ceil(min_value/50);
            for ii=20:-1:complexity_index
                complexity_RL(ii,1)=complexity_RL(ii,1)+1;
            end
        end
    end

end



%探究Q表格的收敛性，这种情况下，每0.5s更新一个距离和连接矩阵，通信拓扑是基于上面的基础上的之后的运动状态。
%仿真跨度为500-3000秒
simulation_time_converge=2500;
time_step_converge=0.5;
% 初始化无人机通信邻接矩阵
comm_matrix_converge = comm_matrix;

%初始化每个离散时间点的距离矩阵
D_converge=cell(simulation_time_converge/time_step_converge,1);
for i=1:(simulation_time_converge/time_step_converge)
    D_converge{i}=zeros(num_drones,num_drones);
end

%初始化每个离散时间点的连接矩阵
I_converge=cell(simulation_time_converge/time_step_converge,1);
for i=1:(simulation_time_converge/time_step_converge)
    I_converge{i}=zeros(num_drones,num_drones);
end

% 模拟无人机运动
for k=1:(simulation_time_converge/time_step_converge)

    %每隔一秒更新无人机速度
    if mod(k,1/time_step_converge)==0
        speed_change=speedChange*randn(num_drones,3);%无人机速度扰动
        drone_speeds=drone_speeds+speed_change;%更新扰动后的速度
    end

    % 更新无人机位置
    drone_positions = drone_positions + drone_speeds * time_step_converge;

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
    D_converge{k}=distances;

    % 更新通信邻接矩阵
    comm_matrix_converge = (distances <= comm_range);
    comm_matrix_converge = comm_matrix_converge-eye(num_drones,num_drones);
    I_converge{k}=comm_matrix_converge;
end

while  min(Tsend(:, 1))<3000
    %先取广播时间最小节点广播，广播完成后按照2,3列数据更新下次广播全局时间
    [min_value, min_index] = min(Tsend(:, 1));
    %获取广播时刻的连接矩阵与距离矩阵
    A_RL=I_converge{ceil((min_value-500)/time_step_converge)};
    B1_RL=D_converge{ceil((min_value-500)/time_step_converge)};

    %更新邻居的频率同步矩阵和相位同步矩阵
    %频率调整第一列为flag位，1表示最新消息，每次同步后要把为1的flag置零。
    %检查min_index节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
    for k=1:num_drones
        if A_RL(min_index,k)==1

            F_RL=Neb_list_fre_RL{k};%频率矩阵
            P_RL=Neb_list_pha_RL{k};%相位矩阵

            %首次接收到的消息（不管是第几轮），存放在前n行，n表示无人机总数
            if F_RL(min_index,2)==0
                F_RL(min_index,1)=1;
                P_RL(min_index,1)=1;
                %填充邻居列表中首次接收到该节点信息的数据(频率和相位矩阵都有)
                for k1=1:5
                    F_RL(min_index,2*k1)=min_value*alpha(min_index,1)+beta(min_index,1)+dt*(k1-1)+5e-9*randn;
                    P_RL(min_index,2*k1)=F_RL(min_index,2*k1);
                end

                for k2=1:5
                    F_RL(min_index,2*k2+1)=(((min_value*alpha(min_index,1)+beta(min_index,1)+(k2-1)*dt-...
                        beta(min_index,1))/(alpha(min_index,1)))+B1_RL(min_index,k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;
                    P_RL(min_index,2*k2+1)=F_RL(min_index,2*k2+1);
                end
                F_RL(min_index,12)=l(min_index,1);
                P_RL(min_index,12)=l(min_index,1);
                P_RL(min_index,13)=h(min_index,1);

                %首次收到相位消息，利用双向消息交换计算距离时延
                T1=F_RL(min_index,2);
                T2=F_RL(min_index,3);
                T3=T2+20e-3;
                T4=(((T3-beta(k,1))/(alpha(k,1)))+B1_RL(min_index,k)/3e8)*alpha(min_index,1)+beta(min_index,1)+5e-9*randn;
                Td=(T2+T4-T3-T1)/2;
                P_RL(min_index,14)=Td;

                Neb_list_fre_RL{k}=F_RL;
                Neb_list_pha_RL{k}=P_RL;
            end

            %非首次最新接受到的消息（不管是那一轮），
            %频率消息存放在第n+1到2n行，用于计算的数据是首次数据和最新数据
            %相位消息首先经过Td补偿，然后更新第一轮
            if F_RL(min_index,2)~=0
                F_RL(min_index+num_drones,1)=1;
                P_RL(min_index,1)=1;
                for k1=1:5
                    F_RL(min_index+num_drones,2*k1)=min_value*alpha(min_index,1)+beta(min_index,1)+dt*(k1-1)+5e-9*randn;
                end

                for k2=1:5
                    F_RL(min_index+num_drones,2*k2+1)=(((min_value*alpha(min_index,1)+beta(min_index,1)+(k2-1)*dt-...
                        beta(min_index,1))/(alpha(min_index,1)))+B1_RL(min_index,k)/3e8)*alpha(k,1)+beta(k,1)+5e-9*randn;
                end
                deltaTd_RL=compensatePha_RL(F_RL,P_RL,min_index,num_drones);%对新消息的相位消息的d进行动态更新
                %更新相位矩阵的值
                for k1=1:5
                    P_RL(min_index,2*k1)=F_RL(min_index+num_drones,2*k1);
                    P_RL(min_index,2*k1+1)=F_RL(min_index+num_drones,2*k1+1);
                end
                P_RL(min_index,12)=l(min_index,1);
                P_RL(min_index,13)=h(min_index,1);
                P_RL(min_index,14)=P_RL(min_index,14)+deltaTd_RL;
                F_RL(min_index+num_drones,12)=l(min_index,1);

                Neb_list_fre_RL{k}=F_RL;
                Neb_list_pha_RL{k}=P_RL;
            end
        end
    end

    %更新Tchange矩阵,用来追踪每个节点每次的动作
    for k3=1:200
        if (Tchange(min_index,k3)==100)
            Tchange(min_index,k3)=Tsend(min_index,3);
            break;
        end
    end

    %第一步：计算方差，获取当前状态，获取上一个动作的奖励，更新上一个状态的Qtable
    P_var_RL=Neb_list_pha_RL{min_index};
    if any(P_var_RL(:) ~= 0) && get_var(P_var_RL,num_drones,l(min_index,1),h(min_index,1))~=0
        var_curr_RL=get_var(P_var_RL,num_drones,l(min_index,1),h(min_index,1));%从相位同步矩阵获取前段时间的时钟值方差，作为当前状态
    else
        var_curr_RL=Var_change(min_index,1);
    end
    var_old_RL=Var_change(min_index,1);%上一轮的方差
    state_old_RL=State_old(min_index,1);%上一轮的转态索引
    %state_curr_RL=0;%这一轮的状态索引初始值

    %获取当前状态索引
    if var_curr_RL-var_old_RL>5e-11 %方差变大
        state_curr_RL=state_to_index(var_curr_RL,1);%状态转索引函数
    elseif var_old_RL-var_curr_RL>5e-11 %方差变小
        state_curr_RL=state_to_index(var_curr_RL,2);
    else %方差几乎保持不变
        state_curr_RL=state_to_index(var_curr_RL,3);
    end

    %将上一轮的方差和状态索引更新为此次的，方便下次使用
    Var_change(min_index,1)=var_curr_RL;
    State_old(min_index,1)=state_curr_RL;

    %从精度，周期，趋势三个方面计算上一个状态采取的动作的奖励R
    R_RL=get_r(Tsend(min_index,3),var_curr_RL,state_curr_RL);

    %Qlearning中的TD方式更新前一个状态的QTable值
    action_old_RL=action_to_index(Tsend(min_index,3));%上一个动作转索引
    max_Q=max(Qtable(state_curr_RL, :));
    Qtable(state_old_RL,action_old_RL)=Qtable(state_old_RL,action_old_RL)+...
        alpha_RL*(R_RL+gamma_RL*max_Q-Qtable(state_old_RL,action_old_RL));

    %第二步：根据贪婪策略选取下次动作，并且更新Tsend的动作，周期和下次发射时间
    action_curr_RL=greedy(state_curr_RL,Qtable);
    Tsend(min_index,3)=action_curr_RL;
    %按照Tsend3列的值更新第一列下次的广播时间和周期
    if Tsend(min_index,2)+Tsend(min_index,3)>2 %周期缩短没有越界
        Tsend(min_index,2)=Tsend(min_index,2)+Tsend(min_index,3);
        Tsend(min_index,1)=(Tsend(min_index,1)*alpha(min_index,1)+beta(min_index,1)+...
            Tsend(min_index,2)-beta(min_index,1))/alpha(min_index,1);
    else
        Tsend(min_index,2)=2;
        Tsend(min_index,1)=(Tsend(min_index,1)*alpha(min_index,1)+beta(min_index,1)+...
            Tsend(min_index,2)-beta(min_index,1))/alpha(min_index,1);
    end

    Qtable_trance{Q_count}=Qtable;
    Q_count=Q_count+1;

    %第三步：完成广播之后，利用自身的频率和相位调整矩阵进行同步
    F_self_RL=Neb_list_fre_RL{min_index};
    P_self_RL=Neb_list_pha_RL{min_index};
    if any(F_self_RL(:) ~= 0) %频率矩阵非空才考虑进行同步
        [l_change,F_change]=distributed_all_l_RL(F_self_RL,num_drones,l(min_index,1));%更新逻辑时钟频偏参数
        %第一个参数表示频偏调整值，第二个表示同步之后的邻居矩阵（要把已经用过的信息标记）
        l(min_index,1)=(1-0.05-0.985^(fix(min_value/2)))*l(min_index,1)+(0.05+0.985^(fix(min_value/2)))*l_change;
        Neb_list_fre_RL{min_index}=F_change;

        [h_change,P_change]=distributed_all_h_RL(P_self_RL,l(min_index,1),h(min_index,1),num_drones);%更新逻辑时钟相偏参数
        h(min_index,1)=(1-0.05-0.985^(fix(min_value/2)))*h(min_index,1)+(0.05+0.985^(fix(min_value/2)))*h_change;%更新逻辑时钟相偏参数
        Neb_list_pha_RL{min_index}=P_change;
    end

    % %这个变量其实看的是400s左右网络的最大时钟差
    % if min(Tsend(:, 1))>390&&min(Tsend(:, 1))<410
    %     global_temp=min(Tsend(:, 1));
    %     local_temp=global_temp*alpha+beta;
    %     compa_temp=local_temp.*l+h;
    %     result_temp=max(compa_temp) - min(compa_temp);
    %     result_zuihou(r,1)=result_temp;
    % end

    %计算强化学习通信开销
        if min_value<=1000
            complexity_index=ceil(min_value/50);
            for ii=20:-1:complexity_index
                complexity_RL(ii,1)=complexity_RL(ii,1)+1;
            end
        end

end

[~, result_converge] = max(Qtable, [], 2);


%五次仿真求平均值
result_test=mean(result_zuihou);
skewFinal=mean(skew);
offsetFinal=mean(offset);
clockFinal=mean(clock);
clock_varFinal=mean(clock_var);
clock_max_eachround(1:25)=clockFinal(1:25);

result_converge_all(3,:)=result_converge;

format long;
%plot(1:simulation_k,skewFinal(1:end),'-*','LineWidth',1);
% hold on;
%plot(1:simulation_k,offsetFinal(1:end),'-*','LineWidth',1);
% hold on;
%semilogy(1:5:simulation_k,clockFinal(1:5:end),'-<','LineWidth',1.5);
%hold on;
non_zero_idx = find(clock_max_eachround ~= 0);  % 获取非零元素的行索引（X轴数据）
non_zero_vals = clock_max_eachround(non_zero_idx);  % 获取对应的非零值（Y轴数据）
is_multiple5 = mod(non_zero_idx, 5) == 0;  % 判断索引是否为5的倍数（逻辑数组）
final_idx = non_zero_idx(is_multiple5);    % 最终X轴：5的倍数的非零索引
final_vals = non_zero_vals(is_multiple5);  % 最终Y轴：对应的值
%semilogy(final_idx,final_vals,'-s','LineWidth',1.5)
% semilogy(1:simulation_k,clock_varFinal(1:end),'-*','LineWidth',1);
%通信开销对比
plot(50:50:1000,complexity_RL(1:20),'-*','LineWidth',1);
hold on;
plot(50:50:1000,complexity_DPTS(1:20),'-s','LineWidth',1);
xlabel('同步轮次');
ylabel('全网最大时钟差/s');
title('密集网络同步性能')
grid on;
hold on;

