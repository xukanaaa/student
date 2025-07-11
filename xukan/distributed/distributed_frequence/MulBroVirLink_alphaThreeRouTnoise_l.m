function M=MulBroVirLink_alphaThreeRouTnoise_l(A,j,i,L)
%采用三轮消息进行偏斜估计，由于此种情况时延（非传播）较大，所以不采用距离补偿算法
% A:节点自身维护的消息矩阵
% j:当前同步轮数
% i:无人机总数
% L:自身节点的频偏补偿值
% l:返回值其实是一个n+1行，2列的矩阵。（1,1）表示同步后的逻辑频偏l,(1,2)表示一共使用的同步数据量
%下面的n行2列矩阵则表示该节点的邻居信息，列名为轮数，斜率αij*lj

%每个节点对应的频率调整值,前i行表示自己的邻居数据，后i*i行则表示每个邻居可能携带的二跳邻居信息
deltal=zeros(i+i*i,1);
%m表示用于同步的总数据量
m=0;
%初始化返回值矩阵
l=zeros(i+1,2);

%遍历矩阵A的下n行
for n=(i+1):(2*i)
    
    %有第三轮该邻居最新消息
    if A(n+i,1)>=(j-1)
        %         %对下n行消息进行偏移d的处理，让其尽量与第一次数据在同一直线上
        %         if A(n,1)-A(n-i,1)<=0 %当两次消息间隔轮次过大时，不需要进行偏移补偿
        %           old_x=(A(n-i,3)+A(n-i,5)+A(n-i,7)+A(n-i,9)+A(n-i,11))/5;%首次接收消息的x轴平均值
        %           old_y=(A(n-i,2)+A(n-i,4)+A(n-i,6)+A(n-i,8)+A(n-i,10))/5;%首次接收消息的y轴平均值
        %           new_x=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11))/5;%最新接收消息的x轴平均值
        %           new_y=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10))/5;%最新接收消息的y轴平均值
        %
        %         %d调整的式子的分子
        %           d_molecule=0;
        %           for k1=1:5
        %               d_molecule=d_molecule+(A(n,2*k1)-old_y)*(new_x-old_x)*(A(n,2*k1+1)-...
        %                   new_x)-(A(n,2*k1+1)-old_x)*(A(n,2*k1+1)-new_x)*(new_y-old_y)...
        %                 -(A(n-i,2*k1+1)-old_x)*(A(n-i,2*k1+1)-old_x)*(new_y-old_y)-...
        %                 (old_y-A(n-i,2*k1))*(A(n-i,2*k1+1)-old_x)*(new_x-old_x);
        %           end
        %
        %         %d调整式子的分母
        %           d_denominator=(A(n,3)-new_x)^2+(A(n,5)-new_x)^2+(A(n,7)-new_x)^2+...
        %               (A(n,9)-new_x)^2+(A(n,11)-new_x)^2+...
        %               (A(n-i,3)-old_x)^2+(A(n-i,5)-old_x)^2+(A(n-i,7)-old_x)^2+...
        %               (A(n-i,9)-old_x)^2+(A(n-i,11)-old_x)^2;
        %
        %         %利用d对最新的邻居消息进行调整（平移）
        %           d=d_molecule/d_denominator;
        %           for k2=1:5
        %               A(n,2*k2)=A(n,2*k2)+d;
        %           end
        %         end
        
        %线性回归求斜率
        Tj_mean=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10)+A(n-i,2)+A(n-i,4)+...
            A(n-i,6)+A(n-i,8)+A(n-i,10)+A(n+i,2)+A(n+i,4)+A(n+i,6)+A(n+i,8)+A(n+i,10))/15;
        Ti_mean=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11)+A(n-i,3)+A(n-i,5)+...
            A(n-i,7)+A(n-i,9)+A(n-i,11)+A(n+i,3)+A(n+i,5)+A(n+i,7)+A(n+i,9)+A(n+i,11))/15;
        
        molecule=0;
        for k3=1:5
            molecule=molecule+(A(n,2*k3)-Tj_mean)*(A(n,2*k3+1)-Ti_mean)+...
                (A(n-i,2*k3)-Tj_mean)*(A(n-i,2*k3+1)-Ti_mean)+(A(n+i,2*k3)-Tj_mean)*(A(n+i,2*k3+1)-Ti_mean);
        end
        denominator=0;
        for k4=1:5
            denominator=denominator+(A(n,2*k4+1)-Ti_mean)^2+(A(n-i,2*k4+1)-Ti_mean)^2+(A(n+i,2*k4+1)-Ti_mean)^2;
        end
        
        deltal(n-i,1)=molecule*A(n+i,12)/denominator;%非首次接收的邻居的斜率，放入deltal中
        %将这个邻居信息加入返回值的邻居列表中
        l(n-i+1,1)=A(n+i,1);
        l(n-i+1,2)=molecule*A(n+i,12)/denominator;
        
        m=m+1;
        
        %遍历该有效邻居信息的二跳矩阵（后n列），计算结果加入deltal中，递增计数器m
        for k5=14:13+i
            if A(n+i,k5)~=0
                deltal(i+i*(n-i-1)+k5-13,1)=(molecule/denominator)*A(n+i,k5);
                m=m+1;
            end
        end
    end
    
    %第二次接收到该邻居消息且是最新消息
    if A(n+i,1)==0 && A(n,1)>=j-1
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
        %将这个邻居信息加入返回值的邻居列表中
        l(n-i+1,1)=A(n,1);
        l(n-i+1,2)=molecule*A(n,12)/denominator;
        
        m=m+1;
        
        %遍历该有效邻居信息的二跳矩阵（后n列），计算结果加入deltal中，递增计数器m
        for k5=14:13+i
            if A(n,k5)~=0
                deltal(i+i*(n-i-1)+k5-13,1)=(molecule/denominator)*A(n,k5);
                m=m+1;
            end
        end
    end
    
    %第一次接收到该邻居消息且是最新消息，不做处理
    %if A(n,1)==0 && A(n-i,1)>=j-1
    %         %利用第一次的消息进行线性回归求斜率
    %        Tj_mean_old=(A(n-i,2)+A(n-i,4)+A(n-i,6)+A(n-i,8)+A(n-i,10))/5;
    %        Ti_mean_old=(A(n-i,3)+A(n-i,5)+A(n-i,7)+A(n-i,9)+A(n-i,11))/5;
    %         molecule=(A(n-i,2)-Tj_mean_old)*(A(n-i,3)-Ti_mean_old)+(A(n-i,4)-Tj_mean_old)*(A(n-i,5)-Ti_mean_old)...
    %             +(A(n-i,6)-Tj_mean_old)*(A(n-i,7)-Ti_mean_old)+(A(n-i,8)-Tj_mean_old)*(A(n-i,9)-Ti_mean_old)+...
    %             (A(n-i,10)-Tj_mean_old)*(A(n-i,11)-Ti_mean_old);
    %         denominator=(A(n-i,3)-Ti_mean_old)^2+(A(n-i,5)-Ti_mean_old)^2+(A(n-i,7)-Ti_mean_old)^2+...
    %             (A(n-i,9)-Ti_mean_old)^2+(A(n-i,11)-Ti_mean_old)^2;
    %
    %         deltal(n-i,1)=molecule*A(n-i,12)/denominator;%把首次接收到消息的节点的更新值加入deltal
    %
    %         %将这个邻居信息加入返回值的邻居列表中
    %         l(n-i+1,1)=A(n-i,1);
    %         l(n-i+1,2)=molecule*A(n-i,12)/denominator;
    %         m=m+1;
    %
    %        %遍历该有效邻居信息的二跳矩阵（后n列），计算结果加入deltal中，递增计数器m
    %        for k5=14:13+i
    %            if A(n-i,k5)~=0
    %                deltal(i+i*(n-i-1)+k5-13,1)=(molecule/denominator)*A(n-i,k5);
    %                m=m+1;
    %            end
    %        end
    %
    %end
end
l(1,1)=(sum(deltal)+L)/(m+1);
l(1,2)=m;
M=l;
end