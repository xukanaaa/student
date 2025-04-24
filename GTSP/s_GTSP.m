num_drones=75; % 无人机数量
space_size =3000; % 空间大小（长、宽、高均为1km）
comm_range =1000; % 通信范围（m）
t_period=1;%发送消息伪周期
simulation_k=400;%仿真轮数
%不同轮次调整后的频偏和相偏
x=zeros(num_drones,simulation_k);
y=zeros(num_drones,simulation_k);
deltamax=zeros(simulation_k,1);%每轮迭代后全网最大钟差
t=zeros(simulation_k,1);%每轮迭代后用于测评的真实时间，此处选为最后一个发送消息的节点每次发送完的时间
speed_change=zeros(num_drones,3);%无人机初始化速度扰动

% 初始化无人机位置
drone_positions = space_size * rand(num_drones, 3);
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
  D=zeros(num_drones,num_drones);


%初始化每个离散时间点的连接矩阵

    I=zeros(num_drones,num_drones);

    
    % 计算无人机之间的距离
    distances = pdist2(drone_positions, drone_positions);
    D=distances;
    
    % 更新通信邻接矩阵
    comm_matrix = (distances <= comm_range);
    comm_matrix = comm_matrix-eye(num_drones,num_drones);
    I=comm_matrix;

t_local_total=zeros(num_drones,simulation_k);
t_global_total=zeros(num_drones,simulation_k);
t_global=(t_period/2)*sort(rand(num_drones,1));%初始发送全球时间，升序
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
for i=1:num_drones
    Neb_list{i}=zeros(num_drones,8);
end

%消息交换开始
for j=1:simulation_k
    for i=1:num_drones
        [c,d]=sort(t_global_total(:,j));
        
        %第一轮不更新，从第二轮开始
        if j>=2
            l(d(i),1)=update_l(Neb_list{d(i)},j,num_drones,l(d(i),1));%更新逻辑时钟频偏参数
            h(d(i),1)=update_h(Neb_list{d(i)},j,num_drones,h(d(i),1));%更新逻辑时钟相偏参数
        end
        
        A=I;%发送消息时刻的连接矩阵
        B=D;%发送消息时刻的距离矩阵
        
        %检查d(i)节点发送时与其有连接的其余节点，更新连接节点的邻居信息列表
        for k=1:num_drones
            if A(d(i),k)==1
                C=Neb_list{k};
                C(d(i),6)=C(d(i),2);
                C(d(i),7)=C(d(i),3);
                C(d(i),1)=j;
                C(d(i),2)=t_local_total(d(i),j)+2e-8*randn;
                C(d(i),3)=(((t_local_total(d(i),j)-beta(d(i),1))/(alpha(d(i),1)))+B(d(i),k)/3e8)*alpha(k,1)+beta(k,1)+2e-8*randn;
                C(d(i),4)=l(d(i),1)*C(d(i),2)+h(d(i),1);
                C(d(i),5)=l(k,1)*C(d(i),3)+h(k,1);
                C(d(i),8)=l(d(i),1);
                Neb_list{k}=C;
            end
        end
    end
  x(:,j)=l.*alpha;
  y(:,j)=l.*beta+h;
end

for i=1:simulation_k
    t(i,1)=(t_local_total(num_drones,i)-beta(num_drones))/alpha(num_drones);
end
for k=1:simulation_k
    deltamax(k,1)=max(x(:,k)*t(k,1)+y(:,k))-min(x(:,k)*t(k,1)+y(:,k));
end
format long;
disp(deltamax)