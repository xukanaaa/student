%根据方差判断网络同步状态，从而动态选取同步邻居数量
function N=distributed_all_vary_h(A,j,num,count,di,l,h)
% A:节点自身维护的消息矩阵
% j:当前同步轮数
% num:无人机总数
% count:频率同步利用的邻居数，相位同步与其一致
% di:当前同步节点编号
% l:频率补偿值矩阵
% h:相位补偿值矩阵
% N:返回值，表示h补偿值

%初始化返回值矩阵
N=0;

%逻辑时钟相位调整
deltah=zeros(num,1);%每个节点对应的相位调整值

m=1;
B=zeros(2,5);%第一行表示邻居的逻辑时钟，第二行表示自身的逻辑时钟

%遍历矩阵A
for n=1:num
    if A(n,1)>=j-1
        for k1=1:5
            A(n,2*k1+1)=A(n,2*k1+1)-A(n,14);%传播时延补偿
        end
        for k1=1:5
            B(1,k1)=A(n,2*k1)*A(n,12)+A(n,13);
            B(2,k1)=A(n,2*k1+1)*l(di)+h(di);
        end
        deltah(n,1)=(B(1,1)-B(2,1)+B(1,2)-B(2,2)+B(1,3)-B(2,3)+B(1,4)-B(2,4)+B(1,5)-B(2,5))/5;
        m=m+1;        
    end
end

%计算与所有邻居时钟差的方差
var_clock=zeros(m-1,5);%初始化方差矩阵
k2=1;
for n=1:num
    if A(n,1)>=j-1
        for n1=1:5
        var_clock(k2,n1)=A(n,2*n1)*A(n,12)+A(n,13)-(A(n,2*n1+1)*l(di)+h(di));
        end
        k2=k2+1;
    end
end
%最终方差
var_result=var(var_clock(:));%最终的方差值
%方差较小，在上一次基础上继续减少同步个数，较大不处理（默认按照所有邻居同步）
if var_result<10e-12
    %减少同步量后仍然大于本次邻居数量，依旧采用当前数量同步
    %减少同步量后小于本次邻居数量，采用减少后的数量同步，具体做法是把deltal中随机m-count+1个非零值置零
    if count<m-1&&count>=1
        % 确定需要置零的元素数量
        k1 = m-count-1;
        
        % 找到所有非零元素的索引
        non_zero_indices = find(deltah ~= 0);
        num_non_zero = length(non_zero_indices);
        
        % 验证非零元素数量是否足够
        if num_non_zero > k1
            % 随机选择k个非零元素的索引
            indices_to_zero = randsample(non_zero_indices, k1);
            
            % 将选中的元素置零
            deltah(indices_to_zero) = 0;
            m=count+1;
        end
    end
end
N=h(di)+(sum(deltah))/m;
end