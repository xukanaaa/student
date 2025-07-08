%使用首轮和最新轮次消息用来同步，并且对相对运动产生的位置偏移进行补偿，利用线性回归最小二乘滤波
%仿真频率同步，相位同步，时钟值同步的性能
%采用虚拟链接的方法加快收敛速度，同时应用于频率同步和相位同步
%同步完成后频偏发生变化，导致同步丧失，但是同步算法不再运行

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
round=2;%蒙特卡洛仿真轮数
skew=zeros(round,simulation_k);%未平均之前的最大频偏差值矩阵
skewFinal=zeros(1,simulation_k);%平均后的每一轮频偏差值。
phySkew=80e-6;%物理时钟频偏最大值
phyOffset=1;%物理时钟相位偏移最大值
offset=zeros(round,simulation_k);
offsetFinal=zeros(1,simulation_k);
clock=zeros(round,simulation_k);
clockFinal=zeros(1,simulation_k);
used_num=0;%每轮同步使用的数据量
M=zeros(num_drones+1,2);

for r=1:2
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
    %不稳定后的物理时钟频偏
    alpha_vary=alpha+(5*randn(num_drones,1)+20)*10^(-6)*((-1)^randi([0,1]));
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
    
    %所有轮次中各节点发送消息的本地物理时间(分两部分，初始偏移表示的部分，和偏移突变之后的部分)
    for m=2:100
        t_local_total(:,m)=t_local+t_period*(m-1)*alpha;
    end
    for m=101:simulation_k
        t_local_total(:,m)=t_local_total(:,100)+t_period*(m-100)*alpha_vary;
    end
    
    %所有轮次中各节点发送消息的全局时间
    for m=2:100
        t_global_total(:,m)=(t_local_total(:,m)-beta)./alpha;
    end
    for m=101:simulation_k
        t_global_total(:,m)=(t_local_total(:,m)-beta)./alpha_vary;
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
    
    for j=1:25
        for i=1:num_drones
            [~,d]=sort(t_global_total(:,j));
                        
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
            
            
                A=I{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的连接矩阵
                B1=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的距离矩阵
                B2=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))+1};        
            
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
    
    for i=26:100
        x(:,i)=l.*alpha;
    end
    for i=101:simulation_k
        x(:,i)=l.*alpha_vary;
    end
    %选取每次评估同步效果的时间（每轮结束评估一次）
    for i=1:100
        t(i,1)=(t_local_total(num_drones,i)-beta(num_drones))/alpha(num_drones);
    end
    for i=101:simulation_k
        t(i,1)=(t_local_total(num_drones,i)-beta(num_drones))/alpha_vary(num_drones);
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
    for d=1:simulation_k
        y1(d,1)=max(y(:,d))-min(y(:,d));
    end
    %一共进行5次蒙特卡洛仿真，r=5，每一行表示各轮次同步后的全网最大频率差
    skew(r,:)=x1';
    offset(r,:)=y1';
    clock(r,:)=deltamax';
end
%五次仿真求平均值
skewFinal=mean(skew);
offsetFinal=mean(offset);
clockFinal=mean(clock);


format long;
% plot(1:simulation_k,skewFinal(1:end),'-*','LineWidth',1);
% hold on;
%plot(1:simulation_k,offsetFinal(1:end),'-*','LineWidth',1);
% hold on;
semilogy(1:simulation_k,clockFinal(1:end),'-mo','LineWidth',1);
xlabel('同步轮数');  % 添加x轴标签
ylabel('网络最大时钟差/秒');      % 添加y轴标签
grid on;
hold on;

