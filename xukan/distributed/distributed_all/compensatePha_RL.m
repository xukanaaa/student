%强化学习自适应周期同步过程中的逻辑时钟频相位调整,对多轮消息的距离进行补偿 
function deltaTd=compensatePha_RL(C,P,k,n)
%对相位消息的距离偏移进行补偿
%C中的第二轮消息表示最新的物理时钟时间戳对
%P中的消息表示上一轮的物理时钟时间戳对
%k表示当前收到节点序号为k的节点的消息
%n表示无人机数目

old_x=(P(k,3)+P(k,5)+P(k,7)+P(k,9)+P(k,11))/5;%首次接收消息的x轴平均值
old_y=(P(k,2)+P(k,4)+P(k,6)+P(k,8)+P(k,10))/5;%首次接收消息的y轴平均值
new_x=(C(n+k,3)+C(n+k,5)+C(n+k,7)+C(n+k,9)+C(n+k,11))/5;%最新接收消息的x轴平均值
new_y=(C(n+k,2)+C(n+k,4)+C(n+k,6)+C(n+k,8)+C(n+k,10))/5;%最新接收消息的y轴平均值

%d调整的式子的分子
d_molecule=0;
for k1=1:5
    d_molecule=d_molecule+(C(n+k,2*k1)-old_y)*(new_x-old_x)*(C(n+k,2*k1+1)-new_x)-...
        (C(n+k,2*k1+1)-old_x)*(C(n+k,2*k1+1)-new_x)*(new_y-old_y)...
        -(P(k,2*k1+1)-old_x)*(P(k,2*k1+1)-old_x)*(new_y-old_y)-(old_y-P(k,2*k1))*(P(k,2*k1+1)-old_x)*(new_x-old_x);
end

%d调整式子的分母
d_denominator=(C(n+k,3)-new_x)^2+(C(n+k,5)-new_x)^2+(C(n+k,7)-new_x)^2+(C(n+k,9)-new_x)^2+(C(n+k,11)-new_x)^2+...
    (P(k,3)-old_x)^2+(P(k,5)-old_x)^2+(P(k,7)-old_x)^2+(P(k,9)-old_x)^2+(P(k,11)-old_x)^2;

%利用d对最新的邻居消息进行调整（平移）
deltaTd=d_molecule/d_denominator;
end