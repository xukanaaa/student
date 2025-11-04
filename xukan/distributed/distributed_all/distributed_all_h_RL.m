%强化学习自适应周期同步过程中的逻辑时钟频相位调整
function [h,P]=distributed_all_h_RL(A,L,H,num)
% A:节点自身维护的消息矩阵
% num:无人机总数
% L:自身节点频率补偿值
% H:自身节点相位补偿值
% h补偿值
% P表示改变后的A矩阵，用于标记已更新数据


%逻辑时钟相位调整
deltah=zeros(num,1);%每个节点对应的相位调整值

m=1;
B=zeros(2,5);%第一行表示邻居的逻辑时钟，第二行表示自身的逻辑时钟
temp_A=A;

%遍历矩阵A
for n=1:num
    if A(n,1)==1
        for k1=1:5
            A(n,2*k1+1)=A(n,2*k1+1)-A(n,14);%传播时延补偿
        end
        for k1=1:5
            B(1,k1)=A(n,2*k1)*A(n,12)+A(n,13);
            B(2,k1)=A(n,2*k1+1)*L+H;
        end
        deltah(n,1)=(B(1,1)-B(2,1)+B(1,2)-B(2,2)+B(1,3)-B(2,3)+B(1,4)-B(2,4)+B(1,5)-B(2,5))/5;
        m=m+1;
        temp_A(n,1)=0;%把使用过的数据标记为旧数据
    end
end
h=H+(sum(deltah))/m;
P=temp_A;
end