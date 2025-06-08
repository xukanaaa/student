function h=MulBro_TwoRouTnoise_h(A,j,L,H,num)
% A:节点自身维护的消息矩阵
% j:当前同步轮数
% num:无人机总数
% L:自身节点频率补偿值
% H:自身节点相位补偿值

%逻辑时钟相位调整
deltah=zeros(num,1);%每个节点对应的相位调整值
% alpha=zeros(5,1);
m=1;
B=zeros(2,5);%第一列表示邻居的逻辑时钟，第二列表示自身的逻辑时钟

%遍历矩阵A的下n行
for n=1:num
    if A(n,1)>=j-1
        for k1=1:5
            A(n,2*k1+1)=A(n,2*k1+1)-A(n,14);%传播时延补偿
        end
        for k1=1:5
            B(1,k1)=A(n,2*k1)*A(n,12)+A(n,13);
            B(2,k1)=A(n,2*k1+1)*L+H;
        end
        deltah(n,1)=(B(1,1)-B(2,1)+B(1,2)-B(2,2)+B(1,3)-B(2,3)+B(1,4)-B(2,4)+B(1,5)-B(2,5))/5;
        m=m+1;
    end
end
h=H+(sum(deltah))/m;
end