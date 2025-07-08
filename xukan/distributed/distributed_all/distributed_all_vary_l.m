function M=distributed_all_vary_l(A,j,i,L,count,di,l,h)
% A:节点自身维护的消息矩阵
% j:当前同步轮数
% i:无人机总数
% L:自身节点的频偏补偿值
% count:表示上一轮计算所用的邻居数量
% di:表示当前是哪个节点更新参数
% l:表示当前逻辑时钟频率调整值，用于计算逻辑时间
% h:表示当前逻辑时间相位调整值，同于计算逻辑时间
% M:返回值矩阵。第一行分别是更新用的l,更新所用的邻居数量


%每个节点对应的频率调整值,前i行表示自己的邻居数据
deltal=zeros(i,1);
%m表示用于同步的总数据量
m=0;
%初始化返回值矩阵
M=zeros(1,2);

%遍历矩阵A的下n行
for n=(i+1):(2*i)   
    %有该邻居最新消息，并且不是第一次获得该邻居信息
    if A(n,1)>=(j-1)
        
        %对下n行消息进行偏移d的处理，让其尽量与第一次数据在同一直线上
        if A(n,1)-A(n-i,1)<=30 %当两次消息间隔轮次过大时，不需要进行偏移补偿
            old_x=(A(n-i,3)+A(n-i,5)+A(n-i,7)+A(n-i,9)+A(n-i,11))/5;%首次接收消息的x轴平均值
            old_y=(A(n-i,2)+A(n-i,4)+A(n-i,6)+A(n-i,8)+A(n-i,10))/5;%首次接收消息的y轴平均值
            new_x=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11))/5;%最新接收消息的x轴平均值
            new_y=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10))/5;%最新接收消息的y轴平均值
            
            %d调整的式子的分子
            d_molecule=0;
            for k1=1:5
                d_molecule=d_molecule+(A(n,2*k1)-old_y)*(new_x-old_x)*(A(n,2*k1+1)-...
                    new_x)-(A(n,2*k1+1)-old_x)*(A(n,2*k1+1)-new_x)*(new_y-old_y)...
                    -(A(n-i,2*k1+1)-old_x)*(A(n-i,2*k1+1)-old_x)*(new_y-old_y)-...
                    (old_y-A(n-i,2*k1))*(A(n-i,2*k1+1)-old_x)*(new_x-old_x);
            end
            
            %d调整式子的分母
            d_denominator=(A(n,3)-new_x)^2+(A(n,5)-new_x)^2+(A(n,7)-new_x)^2+...
                (A(n,9)-new_x)^2+(A(n,11)-new_x)^2+...
                (A(n-i,3)-old_x)^2+(A(n-i,5)-old_x)^2+(A(n-i,7)-old_x)^2+...
                (A(n-i,9)-old_x)^2+(A(n-i,11)-old_x)^2;
            
            %利用d对最新的邻居消息进行调整（平移）
            d=d_molecule/d_denominator;
            for k2=1:5
                A(n,2*k2)=A(n,2*k2)+d;
            end
        end
        
        %线性回归求斜率
        Tj_mean=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10)+A(n-i,2)+A(n-i,4)+...
            A(n-i,6)+A(n-i,8)+A(n-i,10))/10;
        Ti_mean=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11)+A(n-i,3)+A(n-i,5)+...
            A(n-i,7)+A(n-i,9)+A(n-i,11))/10;
        
        molecule=0;
        for k3=1:5
            molecule=molecule+(A(n,2*k3)-Tj_mean)*(A(n,2*k3+1)-Ti_mean)+...
                (A(n-i,2*k3)-Tj_mean)*(A(n-i,2*k3+1)-Ti_mean);
        end
        denominator=0;
        for k4=1:5
            denominator=denominator+(A(n,2*k4+1)-Ti_mean)^2+(A(n-i,2*k4+1)-Ti_mean)^2;
        end
        
        deltal(n-i,1)=molecule*A(n,12)/denominator;%非首次接收的邻居的斜率，放入deltal中
        
        m=m+1;
    end
    
    %第一次接收到该邻居消息且是最新消息
    if A(n,1)==0 && A(n-i,1)>=j-1
        
        %利用第一次的消息进行线性回归求斜率
        Tj_mean_old=(A(n-i,2)+A(n-i,4)+A(n-i,6)+A(n-i,8)+A(n-i,10))/5;
        Ti_mean_old=(A(n-i,3)+A(n-i,5)+A(n-i,7)+A(n-i,9)+A(n-i,11))/5;
        molecule=(A(n-i,2)-Tj_mean_old)*(A(n-i,3)-Ti_mean_old)+(A(n-i,4)-Tj_mean_old)*(A(n-i,5)-Ti_mean_old)...
            +(A(n-i,6)-Tj_mean_old)*(A(n-i,7)-Ti_mean_old)+(A(n-i,8)-Tj_mean_old)*(A(n-i,9)-Ti_mean_old)+...
            (A(n-i,10)-Tj_mean_old)*(A(n-i,11)-Ti_mean_old);
        denominator=(A(n-i,3)-Ti_mean_old)^2+(A(n-i,5)-Ti_mean_old)^2+(A(n-i,7)-Ti_mean_old)^2+...
            (A(n-i,9)-Ti_mean_old)^2+(A(n-i,11)-Ti_mean_old)^2;
        
        deltal(n-i,1)=molecule*A(n-i,12)/denominator;%把首次接收到消息的节点的更新值加入deltal
        
        m=m+1;
    end
end

%对自身和邻居的时钟值求方差，通过方差来判断用来计算的邻居数量应该多还是少
%计算有多少邻居
m1=0;
for n=1:2*i
    if A(n,1)>=(j-1)
        m1=m1+1;
    end
end
%用于计算方差的数据矩阵
var_clock=zeros(m1,5);

%填充数据
k1=1;
for n=1:2*i
    if A(n,1)>=(j-1)
        if n<=i
            for n1=1:5
                var_clock(k1,n1)=A(n,2*n1)*l(n,1)+h(n,1)-(A(n,2*n1+1)*l(di,1)+h(di,1));
            end
            k1=k1+1;
        end
        if n>i
            for n1=1:5
                var_clock(k1,n1)=A(n,2*n1)*l(n-i,1)+h(n-i,1)-(A(n,2*n1+1)*l(di,1)+h(di,1));
            end
             k1=k1+1;
        end
    end
end
var_result=var(var_clock(:));%最终的方差值

%方差较小的情况可以在上一次基础上再减少同步数量
if var_result<10e-12
    %减少同步量后仍然大于本次邻居数量，依旧采用当前数量同步
    %减少同步量后小于本次邻居数量，采用减少后的数量同步，具体做法是把deltal中随机m-count+1个非零值置零
    if count-1<m&&count-1>=1
        % 确定需要置零的元素数量
        k1 = m-count+1;
        
        % 找到所有非零元素的索引
        non_zero_indices = find(deltal ~= 0);
        num_non_zero = length(non_zero_indices);
        
        % 验证非零元素数量是否足够
        if num_non_zero > k1
            % 随机选择k个非零元素的索引
            indices_to_zero = randsample(non_zero_indices, k1);
            
            % 将选中的元素置零
            deltal(indices_to_zero) = 0;
            m=count-1;
        end
    end
    
    if count-1==0&&count-1<m
        % 确定需要置零的元素数量
        k1 = m-1;
        
        % 找到所有非零元素的索引
        non_zero_indices = find(deltal ~= 0);
        num_non_zero = length(non_zero_indices);
        
        % 验证非零元素数量是否足够
        if num_non_zero > k1
            % 随机选择k个非零元素的索引
            indices_to_zero = randsample(non_zero_indices, k1);
            
            % 将选中的元素置零
            deltal(indices_to_zero) = 0;
            m=1;
        end
    end
    %方差较大的情况，直接利用所有信息同步
end

M(1,1)=(sum(deltal)+L)/(m+1);
M(1,2)=m;
end