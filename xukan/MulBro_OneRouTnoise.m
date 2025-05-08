%只使用一次广播的五条消息进行同步，配合最小二乘线性回归滤波，但是由于时间间隔短，所以噪声对斜率估计的影响大。

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

%初始化每个节点本地的邻居信息列表(消息轮，(当前邻居物理，当前自身物理)*广播轮数，当前邻居逻辑时钟频率调整l)
Neb_list=cell(num_drones,1);
for i=1:num_drones
    Neb_list{i}=zeros(num_drones,broad*2+2);
end

for j=1:simulation_k
    for i=1:num_drones
        [c,d]=sort(t_global_total(:,j));
        
             %第一轮不更新，从第二轮开始
        if j>=2
            l(d(i),1)=MulBro_OneRouTnoise_l(Neb_list{d(i)},j,num_drones,l(d(i),1));%更新逻辑时钟频偏参数
            %h(d(i),1)=update_h(Neb_list{d(i)},j,num_drones,h(d(i),1));%更新逻辑时钟相偏参数
        end
        
        A=I{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的连接矩阵
        B1=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的距离矩阵
        B2=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))+1};
       
        %检查d(i)节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
        for k=1:num_drones
            if A(d(i),k)==1
                C=Neb_list{k};
                C(d(i),1)=j;
                %计算一次同步的五次广播的发送和接收时间
                for k1=1:5
                    C(d(i),2*k1)=t_local_total(d(i),j)+dt*(k1-1)+8e-7*randn;
                end
                
                for k1=1:5
                    C(d(i),2*k1+1)=(((t_local_total(d(i),j)+dt*(k1-1)-beta(d(i),1))/(alpha(d(i),1)))+B1(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+8e-7*randn;
                end
              
                C(d(i),12)=l(d(i),1);
                Neb_list{k}=C;
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
end
%五次仿真求平均值
skewFinal=mean(skew);

format long;
plot(20:simulation_k,skewFinal(20:end),'-*');

