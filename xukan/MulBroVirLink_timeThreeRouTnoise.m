%使用首轮和最新轮次消息用来同步，并且对相对运动产生的位置偏移进行补偿，利用线性回归最小二乘滤波
%采用选择性区域性虚拟链接，增加图的代数联通度，从而加快收敛速度，多跳虚拟链接传递的不是斜率，而是时间戳的值
%采用三轮消息进行估计，增加精确度。

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
t=zeros(simulation_k,1);%每轮迭代后用于测评的真实时间，此处选为最后一个发送消息的节点每次发送完的时间
speed_change=zeros(num_drones,3);%无人机初始化速度扰动
round=5;%蒙特卡洛仿真轮数
skew=zeros(round,simulation_k);%未平均之前的最大频偏差值矩阵
skewFinal=zeros(1,simulation_k);%平均后的每一轮频偏差值。
phySkew=80e-6;%物理时钟频偏最大值
phyOffset=1;%物理时钟相位偏移最大值
used_num=0;%每轮同步使用的数据量
M=zeros(num_drones+1,2);

for r=1:5
    
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
    
    %初始化每个节点本地的用于同步邻居信息列表(消息轮，(当前邻居物理，当前自身物理)*广播轮数，当前邻居逻辑时钟频率调整l,
    %flag位，只占该列第一行，下面都是0,多跳信息是否覆盖位（测量值为0，多跳信息由一轮数据估计为1,2轮为2,3轮为3,0覆盖123,3覆盖12,2覆盖1等等）
    %前n行表示首次接收到某节点的消息数据
    %中n行表示最新或者中间接收到某节点的消息数据
    %后n行表示最新接收到某节点的消息数据
    Neb_list=cell(num_drones,1);
    for i=1:num_drones
        Neb_list{i}=zeros(num_drones*3,1+broad*2+2+1);
    end
    %初始化每个节点本地的每一轮接收到的邻居信息列表(消息轮，(当前邻居物理，当前自身物理)*广播轮数（不包含多跳邻居）
    %用于计算当前节点所处网络的稠密程度
    Neb_list1=cell(num_drones,1);
    for i=1:num_drones
        Neb_list1{i}=zeros(num_drones,1+broad*2);
    end
    
    %模拟同步更新和广播
    for j=1:simulation_k
        for i=1:num_drones
            [c,d]=sort(t_global_total(:,j));
            
            %第一轮不更新，从第二轮开始
            if j>=2
                l(d(i),1)=MulBroVirLink_timeThreeRouTnoise_l(Neb_list{d(i)},j,num_drones,l(d(i),1));%更新逻辑时钟频偏参数
            end
            
            A=I{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的连接矩阵
            B1=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的距离矩阵
            B2=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))+1};
            
            
            D1=Neb_list{d(i)};%广播节点的自身的用于同步的消息矩阵
            D2=Neb_list1{d(i)};%广播节点每轮同步的邻居信息，用于判断局部网络稀疏性
            countNeb=0;
            for k=1:num_drones
                if(D2(k,1)>=j-1)
                    countNeb=countNeb+1;
                end
            end
            %检查d(i)节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
            for k=1:num_drones
                if A(d(i),k)==1
                    C=Neb_list{k};%邻居的消息矩阵,用于同步
                    C1=Neb_list1{k};%邻居的每轮同步信息，用于判断网络稀疏性
                    
                    %对于自己上一轮用于同步邻居小于3时，向所有邻居发送求助，修改邻居的flag为1，
                    %该邻居检测到flag为1时，则在下次广播中附加自己的二跳信息
                    if j>=2 && countNeb<=200
                        C(1,13)=1;
                    end
                    
                    %首次接收到的消息（不管是第几轮）或者第一轮消息为多跳，一律替换成单跳，存放在前n行，n表示无人机总数
                    if C(d(i),1)==0 || C(d(i),14)==1||C(d(i),14)==2||C(d(i),14)==3
                        C(d(i),1)=j;
                        C1(d(i),1)=j;
                        %填充邻居列表中首次接收到该节点信息的数据,更新用于同步的表和单纯邻居表
                        for k1=1:5
                            C(d(i),2*k1)=t_local_total(d(i),j)+dt*(k1-1)+2e-6*randn+2e-8*randn;
                            C1(d(i),2*k1)=C(d(i),2*k1);
                        end
                        
                        for k2=1:5
                            C(d(i),2*k2+1)=(((t_local_total(d(i),j)+(k2-1)*dt-...
                                beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+2e-8*randn;
                            C1(d(i),2*k2+1)=C(d(i),2*k2+1);
                        end
                        C(d(i),12)=l(d(i),1);
                        C(d(i),14)=0;%单跳广播得到的消息最后一个标志位为0
                        
                        %判断是否往邻居消息列表中加入当前广播节点的二跳信息，如果广播节点flag=1，则需要添加
                        %添加规则是，如果（第一轮为空)或者（第一轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                        %如果第二轮为空or第二轮是最新消息也是二跳且估计数据集小于当前二跳数据集，则覆盖
                        %如果第三轮为空or第三轮不是最新消息or第三轮是最新消息也是二跳且二跳数据集小于当前二跳估计数据集，则覆盖
                        %这样，同一轮可能受到同一个二跳节点信息，取其中较准确的值。前两轮只要有数据就不覆盖，只覆盖准确性较低的二跳数据或者空
                        %最后一轮如果没有最新数据就覆盖，考虑一跳，二跳等情况
                        if D1(1,13)==1
                            %遍历邻居列表，找出上一轮同步的邻居
                            for k3=1:num_drones
                                if D2(k3,1)>=j-1
                                    %对上一轮同步的邻居的时钟值进行补偿，返回补偿后的时钟值和所用的消息轮数
                                    deltat=C(d(i),2)-D2(k3,3);%广播节点上次接收到消息到此次广播消息的时间间隔
                                    [alphat,m1]=compensate(D1,k3,num_drones,j);%返回斜率和估计所用的轮数
                                    deltat_y=alphat*deltat;%用于对二跳邻居时钟补偿的时间差
                                    
                                    %k3对应的二跳邻居信息是否可以加入一轮同步矩阵
                                    %（第一轮为空)或者（第一轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                    if (C(k3,1)==0)||(C(k3,1)==j&&C(k3,14)>0&&C(k3,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3,1)=j;
                                        for k4=1:5
                                            C(k3,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3,2*k4+1)=C(d(i),2*k4+1);
                                        end
                                        C(k3,14)=m1;
                                        C(k3,12)=l(k3,1);
                                        %k3对应的二跳邻居信息是否可以加入二轮同步矩阵
                                        %（第二轮为空)或者（第二轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                    elseif (C(k3+num_drones,1)==0)||(C(k3+num_drones,1)==j&&C(k3+num_drones,14)>0&&C(k3+num_drones,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3+num_drones,1)=j;
                                        for k4=1:5
                                            C(k3+num_drones,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3+num_drones,2*k4+1)=C(d(i),2*k4+1);
                                        end
                                        C(k3+num_drones,14)=m1;
                                        C(k3+num_drones,12)=l(k3,1);
                                        %如果第三轮不是最新消息or第三轮是最新消息也是二跳且二跳数据集小于当前二跳估计数据集，则覆盖
                                    elseif (C(k3+num_drones*2,1)<j)||(C(k3+num_drones*2,1)==j&&C(k3+num_drones*2,14)>0&&C(k3+num_drones*2,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3+num_drones*2,1)=j;
                                        for k4=1:5
                                            C(k3+num_drones*2,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3+num_drones*2,2*k4+1)=C(d(i),2*k4+1);
                                        end
                                        C(k3+num_drones*2,14)=m1;
                                        C(k3+num_drones*2,12)=l(k3,1);
                                    end
                                end
                            end
                        end
                        Neb_list{k}=C;
                        Neb_list1{k}=C1;
                        
                        
                        %非首次最新接受到的消息（第二轮收到），存放在第n+1到2n行，用于计算的数据是首次数据和最新数据
                    elseif C(d(i)+num_drones,1)==0||C(d(i)+num_drones,14)==1||C(d(i)+num_drones,14)==2||C(d(i)+num_drones,14)==3
                        C(d(i)+num_drones,1)=j;
                        C1(d(i)+num_drones,1)=j;
                        for k1=1:5
                            C(d(i)+num_drones,2*k1)=t_local_total(d(i),j)+dt*(k1-1)+2e-8*randn+2e-6*randn;
                            C1(d(i)+num_drones,2*k1)=C(d(i)+num_drones,2*k1);
                        end
                        
                        for k2=1:5
                            C(d(i)+num_drones,2*k2+1)=(((t_local_total(d(i),j)+(k2-1)*dt-...
                                beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+2e-8*randn;
                            C1(d(i)+num_drones,2*k2+1)=C(d(i)+num_drones,2*k2+1);
                        end
                        C(d(i)+num_drones,12)=l(d(i),1);
                        C(d(i)+num_drones,14)=0;%单跳广播得到的消息最后一个标志位为0
                        %判断是否往邻居消息列表中加入当前广播节点的二跳信息，如果广播节点flag=1，则需要添加
                        if D1(1,13)==1
                            %遍历邻居列表，找出上一轮同步的邻居
                            for k3=1:num_drones
                                if D2(k3,1)>=j-1
                                    %对上一轮同步的邻居的时钟值进行补偿，返回补偿后的时钟值和所用的消息轮数
                                    deltat=C(d(i)+num_drones,2)-D2(k3,3);%广播节点上次接收到消息到此次广播消息的时间间隔
                                    [alphat,m1]=compensate(D1,k3,num_drones,j);%返回斜率和估计所用的轮数
                                    deltat_y=alphat*deltat;%用于对二跳邻居时钟补偿的时间差
                                    
                                    %k3对应的二跳邻居信息是否可以加入一轮同步矩阵
                                    %（第一轮为空)或者（第一轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                    if (C(k3,1)==0)||(C(k3,1)==j&&C(k3,14)>0&&C(k3,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3,1)=j;
                                        for k4=1:5
                                            C(k3,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3,2*k4+1)=C(d(i)+num_drones,2*k4+1);
                                        end
                                        C(k3,14)=m1;
                                        C(k3,12)=l(k3,1);
                                        %k3对应的二跳邻居信息是否可以加入二轮同步矩阵
                                        %（第二轮为空)或者（第二轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                    elseif (C(k3+num_drones,1)==0)||(C(k3+num_drones,1)==j&&C(k3+num_drones,14)>0&&C(k3+num_drones,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3+num_drones,1)=j;
                                        for k4=1:5
                                            C(k3+num_drones,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3+num_drones,2*k4+1)=C(d(i)+num_drones,2*k4+1);
                                        end
                                        C(k3+num_drones,14)=m1;
                                        C(k3+num_drones,12)=l(k3,1);
                                        %如果第三轮不是最新消息or第三轮是最新消息也是二跳且二跳数据集小于当前二跳估计数据集，则覆盖
                                    elseif (C(k3+num_drones*2,1)<j)||(C(k3+num_drones*2,1)==j&&C(k3+num_drones*2,14)>0&&C(k3+num_drones*2,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3+num_drones*2,1)=j;
                                        for k4=1:5
                                            C(k3+num_drones*2,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3+num_drones*2,2*k4+1)=C(d(i)+num_drones,2*k4+1);
                                        end
                                        C(k3+num_drones*2,14)=m1;
                                        C(k3+num_drones*2,12)=l(k3,1);
                                    end
                                end
                            end
                        end
                        Neb_list{k}=C;
                        Neb_list1{k}=C1;
                        
                        %第三轮接收到消息，存放在第2n+1到3n行，用于计算的数据是首次数据和最新数据和中间数据
                    elseif C(d(i)+num_drones*2,1)==0||C(d(i)+num_drones*2,14)==1||C(d(i)+num_drones*2,14)==2||C(d(i)+num_drones*2,14)==3
                        C(d(i)+num_drones*2,1)=j;
                        C1(d(i)+num_drones*2,1)=j;
                        for k1=1:5
                            C(d(i)+num_drones*2,2*k1)=t_local_total(d(i),j)+dt*(k1-1)+2e-8*randn+2e-6*randn;
                            C1(d(i)+num_drones*2,2*k1)=C(d(i)+num_drones*2,2*k1);
                        end
                        
                        for k2=1:5
                            C(d(i)+num_drones*2,2*k2+1)=(((t_local_total(d(i),j)+(k2-1)*dt-...
                                beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+2e-8*randn;
                            C1(d(i)+num_drones*2,2*k2+1)=C(d(i)+num_drones*2,2*k2+1);
                        end
                        C(d(i)+num_drones*2,12)=l(d(i),1);
                        C(d(i)+num_drones*2,14)=0;%单跳广播得到的消息最后一个标志位为0
                        
                        %判断是否往邻居消息列表中加入当前广播节点的二跳信息，如果广播节点flag=1，则需要添加
                        if D1(1,13)==1
                            %遍历邻居列表，找出上一轮同步的邻居
                            for k3=1:num_drones
                                if D2(k3,1)>=j-1
                                    %对上一轮同步的邻居的时钟值进行补偿，返回补偿后的时钟值和所用的消息轮数
                                    deltat=C(d(i)+num_drones*2,2)-D2(k3,3);%广播节点上次接收到消息到此次广播消息的时间间隔
                                    [alphat,m1]=compensate(D1,k3,num_drones,j);%返回斜率和估计所用的轮数
                                    deltat_y=alphat*deltat;%用于对二跳邻居时钟补偿的时间差
                                    
                                    %k3对应的二跳邻居信息是否可以加入一轮同步矩阵
                                    %（第一轮为空)或者（第一轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                    if (C(k3,1)==0)||(C(k3,1)==j&&C(k3,14)>0&&C(k3,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3,1)=j;
                                        for k4=1:5
                                            C(k3,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                        end
                                        C(k3,14)=m1;
                                        C(k3,12)=l(k3,1);
                                        %k3对应的二跳邻居信息是否可以加入二轮同步矩阵
                                        %（第二轮为空)或者（第二轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                    elseif (C(k3+num_drones,1)==0)||(C(k3+num_drones,1)==j&&C(k3+num_drones,14)>0&&C(k3+num_drones,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3+num_drones,1)=j;
                                        for k4=1:5
                                            C(k3+num_drones,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3+num_drones,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                        end
                                        C(k3+num_drones,14)=m1;
                                        C(k3+num_drones,12)=l(k3,1);
                                        %如果第三轮不是最新消息or第三轮是最新消息也是二跳且二跳数据集小于当前二跳估计数据集，则覆盖
                                    elseif (C(k3+num_drones*2,1)<j)||(C(k3+num_drones*2,1)==j&&C(k3+num_drones*2,14)>0&&C(k3+num_drones*2,14)<m1)
                                        %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                        C(k3+num_drones*2,1)=j;
                                        for k4=1:5
                                            C(k3+num_drones*2,2*k4)=D2(k3,2*k4)+deltat_y;
                                        end
                                        for k4=1:5
                                            C(k3+num_drones*2,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                        end
                                        C(k3+num_drones*2,14)=m1;
                                        C(k3+num_drones*2,12)=l(k3,1);
                                    end
                                end
                            end
                        end
                        
                        Neb_list{k}=C;
                        Neb_list1{k}=C1;
                        
                        %第三轮接收到消息，并且已经存放了三轮消息，需要对消息记录进行清理。用于计算的数据是首次数据和最新数据和中间数据
                    else
                        
                        %如果最新轮次减中间轮次小于中间轮次减首次，直接覆盖第三轮内容
                        if (C(d(i)+num_drones*2,1)-C(d(i)+num_drones,1))<=(C(d(i)+num_drones,1)-C(d(i),1))
                            C(d(i)+num_drones*2,1)=j;
                            C1(d(i)+num_drones*2,1)=j;
                            for k1=1:5
                                C(d(i)+num_drones*2,2*k1)=t_local_total(d(i),j)+dt*(k1-1)+2e-8*randn+2e-6*randn;
                                C1(d(i)+num_drones*2,2*k1)=C(d(i)+num_drones*2,2*k1);
                            end
                            
                            for k2=1:5
                                C(d(i)+num_drones*2,2*k2+1)=(((t_local_total(d(i),j)+(k2-1)*dt-...
                                    beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+2e-8*randn;
                                C1(d(i)+num_drones*2,2*k2+1)=C(d(i)+num_drones*2,2*k2+1);
                            end
                            C(d(i)+num_drones*2,12)=l(d(i),1);
                            C(d(i)+num_drones*2,14)=0;
                            
                            %判断是否往邻居消息列表中加入当前广播节点的二跳信息，如果广播节点flag=1，则需要添加
                            if D1(1,13)==1
                                %遍历邻居列表，找出上一轮同步的邻居
                                for k3=1:num_drones
                                    if D2(k3,1)>=j-1
                                        %对上一轮同步的邻居的时钟值进行补偿，返回补偿后的时钟值和所用的消息轮数
                                        deltat=C(d(i)+num_drones*2,2)-D2(k3,3);%广播节点上次接收到消息到此次广播消息的时间间隔
                                        [alphat,m1]=compensate(D1,k3,num_drones,j);%返回斜率和估计所用的轮数
                                        deltat_y=alphat*deltat;%用于对二跳邻居时钟补偿的时间差
                                        
                                        %k3对应的二跳邻居信息是否可以加入一轮同步矩阵
                                        %（第一轮为空)或者（第一轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                        if (C(k3,1)==0)||(C(k3,1)==j&&C(k3,14)>0&&C(k3,14)<m1)
                                            %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                            C(k3,1)=j;
                                            for k4=1:5
                                                C(k3,2*k4)=D2(k3,2*k4)+deltat_y;
                                            end
                                            for k4=1:5
                                                C(k3,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                            end
                                            C(k3,14)=m1;
                                            C(k3,12)=l(k3,1);
                                            %k3对应的二跳邻居信息是否可以加入二轮同步矩阵
                                            %（第二轮为空)或者（第二轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                        elseif (C(k3+num_drones,1)==0)||(C(k3+num_drones,1)==j&&C(k3+num_drones,14)>0&&C(k3+num_drones,14)<m1)
                                            %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                            C(k3+num_drones,1)=j;
                                            for k4=1:5
                                                C(k3+num_drones,2*k4)=D2(k3,2*k4)+deltat_y;
                                            end
                                            for k4=1:5
                                                C(k3+num_drones,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                            end
                                            C(k3+num_drones,14)=m1;
                                            C(k3+num_drones,12)=l(k3,1);
                                            %如果第三轮不是最新消息or第三轮是最新消息也是二跳且二跳数据集小于当前二跳估计数据集，则覆盖
                                        elseif (C(k3+num_drones*2,1)<j)||(C(k3+num_drones*2,1)==j&&C(k3+num_drones*2,14)>0&&C(k3+num_drones*2,14)<m1)
                                            %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                            C(k3+num_drones*2,1)=j;
                                            for k4=1:5
                                                C(k3+num_drones*2,2*k4)=D2(k3,2*k4)+deltat_y;
                                            end
                                            for k4=1:5
                                                C(k3+num_drones*2,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                            end
                                            C(k3+num_drones*2,14)=m1;
                                            C(k3+num_drones*2,12)=l(k3,1);
                                        end
                                    end
                                end
                            end
                            Neb_list{k}=C;
                            Neb_list1{k}=C1;
                        end
                        
                        %如果最新轮次减中间轮次大于中间轮次减第一轮，则需要把第三轮数据下放到第二轮，最新数据放到第三轮
                        if (C(d(i)+num_drones*2,1)-C(d(i)+num_drones,1))>(C(d(i)+num_drones,1)-C(d(i),1))
                            %第三轮拷贝到第二轮
                            C(d(i)+num_drones,:)=C(d(i)+num_drones*2,:);
                            %最新轮放到第三轮
                            C(d(i)+num_drones*2,1)=j;
                            C1(d(i)+num_drones*2,1)=j;
                            for k1=1:5
                                C(d(i)+num_drones*2,2*k1)=t_local_total(d(i),j)+dt*(k1-1)+2e-8*randn+2e-6*randn;
                                C1(d(i)+num_drones*2,2*k1)=C(d(i)+num_drones*2,2*k1);
                            end
                            
                            for k2=1:5
                                C(d(i)+num_drones*2,2*k2+1)=(((t_local_total(d(i),j)+(k2-1)*dt-...
                                    beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+2e-8*randn;
                                C1(d(i)+num_drones*2,2*k2+1)=C(d(i)+num_drones*2,2*k2+1);
                            end
                            C(d(i)+num_drones*2,12)=l(d(i),1);
                            C(d(i)+num_drones*2,14)=0;
                            %判断是否往邻居消息列表中加入当前广播节点的二跳信息，如果广播节点flag=1，则需要添加
                            if D1(1,13)==1
                                %遍历邻居列表，找出上一轮同步的邻居
                                for k3=1:num_drones
                                    if D2(k3,1)>=j-1
                                        %对上一轮同步的邻居的时钟值进行补偿，返回补偿后的时钟值和所用的消息轮数
                                        deltat=C(d(i)+num_drones*2,2)-D2(k3,3);%广播节点上次接收到消息到此次广播消息的时间间隔
                                        [alphat,m1]=compensate(D1,k3,num_drones,j);%返回斜率和估计所用的轮数
                                        deltat_y=alphat*deltat;%用于对二跳邻居时钟补偿的时间差
                                        
                                        %k3对应的二跳邻居信息是否可以加入一轮同步矩阵
                                        %（第一轮为空)或者（第一轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                        if (C(k3,1)==0)||(C(k3,1)==j&&C(k3,14)>0&&C(k3,14)<m1)
                                            %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                            C(k3,1)=j;
                                            for k4=1:5
                                                C(k3,2*k4)=D2(k3,2*k4)+deltat_y;
                                            end
                                            for k4=1:5
                                                C(k3,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                            end
                                            C(k3,14)=m1;
                                            C(k3,12)=l(k3,1);
                                            %k3对应的二跳邻居信息是否可以加入二轮同步矩阵
                                            %（第二轮为空)或者（第二轮是最新消息也是二跳且估计数据集小于当前二跳估计数据集），则覆盖
                                        elseif (C(k3+num_drones,1)==0)||(C(k3+num_drones,1)==j&&C(k3+num_drones,14)>0&&C(k3+num_drones,14)<m1)
                                            %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                            C(k3+num_drones,1)=j;
                                            for k4=1:5
                                                C(k3+num_drones,2*k4)=D2(k3,2*k4)+deltat_y;
                                            end
                                            for k4=1:5
                                                C(k3+num_drones,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                            end
                                            C(k3+num_drones,14)=m1;
                                            C(k3+num_drones,12)=l(k3,1);
                                            %如果第三轮不是最新消息or第三轮是最新消息也是二跳且二跳数据集小于当前二跳估计数据集，则覆盖
                                        elseif (C(k3+num_drones*2,1)<j)||(C(k3+num_drones*2,1)==j&&C(k3+num_drones*2,14)>0&&C(k3+num_drones*2,14)<m1)
                                            %分别填充同步轮次，二跳邻居和自身时钟，估计所用数据轮数，l值
                                            C(k3+num_drones*2,1)=j;
                                            for k4=1:5
                                                C(k3+num_drones*2,2*k4)=D2(k3,2*k4)+deltat_y;
                                            end
                                            for k4=1:5
                                                C(k3+num_drones*2,2*k4+1)=C(d(i)+num_drones*2,2*k4+1);
                                            end
                                            C(k3+num_drones*2,14)=m1;
                                            C(k3+num_drones*2,12)=l(k3,1);
                                        end
                                    end
                                end
                            end
                            Neb_list{k}=C;
                            Neb_list1{k}=C1;
                        end
                    end
                end
            end
            
            %当所有的消息都广播完毕后，广播节点的flag要重新置零
            D1(1,13)=0;
            Neb_list{d(i)}=D1;
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
    
    x1=zeros(simulation_k,1);
    y1=zeros(simulation_k,1);
    %每轮结束后的频率偏差
    for c=1:simulation_k
        x1(c,1)=max(x(:,c))-min(x(:,c));
    end
    %每轮结束后的相位偏差
    %for d=1:simulation_k
    %    y1(d,1)=max(y(:,d))-min(y(:,d));
    %end
    %一共进行5次蒙特卡洛仿真，r=5，每一行表示各轮次同步后的全网最大频率差
    skew(r,:)=x1';
end
%五次仿真求平均值
skewFinal=mean(skew);

format long;
plot(5:simulation_k,skewFinal(5:end),'-*');

