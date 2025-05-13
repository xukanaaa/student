num_drones = 100; % 无人机数量
space_size = 9000; % 空间大小（长、宽、高均为1km）
min_speed =80; % 最小速度（m/s）
max_speed =100; % 最大速度（m/s）
comm_range = 3500; % 通信范围（m）
time_step = 0.1; % 时间步长（s）
simulation_time = 500; % 模拟时间（s）
t_period=2;%发送消息伪周期
simulation_k=200;%仿真轮数
%不同轮次调整后的频偏和相偏
x=zeros(num_drones,simulation_k);
y=zeros(num_drones,simulation_k);
z=zeros(num_drones,simulation_k);
deltamax=zeros(simulation_k,1);%每轮迭代后全网最大钟差
t=zeros(simulation_k,1);%每轮迭代后用于测评的真实时间，此处选为最后一个发送消息的节点每次发送完的时间
speed_change=zeros(num_drones,3);%无人机初始化速度扰动
round=5;%仿真轮数
skew=zeros(round,simulation_k);%未平均之前的最大频偏差值矩阵
offset=zeros(round,simulation_k);
skewFinal=zeros(1,simulation_k);%平均后的每一轮频偏差值。
offsetFinal=zeros(1,simulation_k);

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
    % 确保速度向量的模在4-8m/s范围内（调整方向但不改变模）
    %drone_speeds = drone_speeds ./ sqrt(sum(drone_speeds.^2, 2)) .* (min_speed + (max_speed - min_speed) * rand(num_drones, 1));
    %初始化物理时钟频偏
    alpha=1-80e-6+(2*80e-6)*rand(num_drones,1);
    %初始化物理时钟相偏
    beta=-1+2*rand(num_drones,1);
    %初始化逻辑时钟频偏
    l=ones(num_drones,1);
    %初始化逻辑时钟相偏
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
        if mod(k,10)==0
            speed_change=2*randn(num_drones,3);%无人机速度扰动
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
    t_global=(t_period*0.9)*sort(rand(num_drones,1));%初始发送全球时间，升序,同时预留0.1倍周期作为保护间隔
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
    
    %初始化每个节点本地的邻居信息列表(消息轮，当前邻居物理，当前自身物理，当前邻居逻辑，当前自身逻辑，之前邻居物理，之前自身物理，当前邻居l)
    Neb_list=cell(num_drones,1);
    for i=1:num_drones
        Neb_list{i}=zeros(num_drones,8);
    end
    
    %消息交换开始
    for j=1:simulation_k
        for i=1:num_drones
            [c,d]=sort(t_global_total(:,j));
            
            %第一轮不更新，从第二轮开始
            if j>=2
                l(d(i),1)=GTSP_update_l(Neb_list{d(i)},j,num_drones,l(d(i),1));%更新逻辑时钟频偏参数
                h(d(i),1)=GTSP_update_h(Neb_list{d(i)},j,num_drones,h(d(i),1));%更新逻辑时钟相偏参数
            end
            
            A=I{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的连接矩阵
            B=D{ceil((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)*time_step))};%发送消息时刻的距离矩阵
            
            %检查d(i)节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
            for k=1:num_drones
                if A(d(i),k)==1
                    C=Neb_list{k};
                    C(d(i),6)=C(d(i),2);
                    C(d(i),7)=C(d(i),3);
                    C(d(i),1)=j;
                    C(d(i),2)=t_local_total(d(i),j)+2e-9*randn;%频率同步
                    C(d(i),3)=(((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)))+B(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+2e-9*randn;
                    C(d(i),4)=l(d(i),1)*C(d(i),2)+h(d(i),1);%相位同步
                    C(d(i),5)=l(k,1)*C(d(i),3)+h(k,1);
                    C(d(i),8)=l(d(i),1);
                    Neb_list{k}=C;
                end
            end
        end
        x(:,j)=l.*alpha;
        y(:,j)=l.*beta+h;
        z(:,j)=h;
    end
    for i=1:simulation_k
        t(i,1)=(t_local_total(num_drones,i)-beta(num_drones))/alpha(num_drones);
    end
    for k=1:simulation_k
        deltamax(k,1)=max(x(:,k)*t(k,1)+y(:,k))-min(x(:,k)*t(k,1)+y(:,k));
    end
    
    x1=zeros(simulation_k,1);
    y1=zeros(simulation_k,1);
    for c=1:simulation_k
        x1(c,1)=max(x(:,c))-min(x(:,c));
    end
    for d=1:simulation_k
        y1(d,1)=max(y(:,d))-min(y(:,d));
    end
    
    skew(r,:)=x1';
    offset(r,:)=y1';
end

skewFinal=mean(skew);
offsetFinal=mean(offset);

format long;
%plot(5:simulation_k,skewFinal(5:end),'-*','LineWidth',1);
plot(12:simulation_k,offsetFinal(12:end),'-*','LineWidth',1);
xlabel('同步轮次');
ylabel('最大相位偏差');
grid on;
hold on;

%disp(x1)