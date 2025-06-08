function l=MulBro_TwoRouTnoise_l(A,j,i,L)
% A:节点自身维护的消息矩阵
% j:当前同步轮数
% i:无人机总数
% L:自身节点的频偏补偿值

%GTSP的逻辑时钟频率调整
deltal=zeros(i,1);%每个节点对应的频率调整值
% alpha=zeros(5,1);
m=0;

%遍历矩阵A的下n行
for n=(i+1):(2*i)
    %有该邻居最新消息，并且不是第一次获得该邻居信息
    if A(n,1)>=(j-1)
        %对下n行消息进行偏移d的处理，让其尽量与第一次数据在同一直线上
        if A(n,1)-A(n-i,1)<=100
          old_x=(A(n-i,3)+A(n-i,5)+A(n-i,7)+A(n-i,9)+A(n-i,11))/5;%首次接收消息的x轴平均值
          old_y=(A(n-i,2)+A(n-i,4)+A(n-i,6)+A(n-i,8)+A(n-i,10))/5;%首次接收消息的y轴平均值
          new_x=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11))/5;%最新接收消息的x轴平均值
          new_y=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10))/5;%最新接收消息的y轴平均值
        
        %d调整的式子的分子
          d_molecule=0;
          for k1=1:5
              d_molecule=d_molecule+(A(n,2*k1)-old_y)*(new_x-old_x)*(A(n,2*k1+1)-new_x)-(A(n,2*k1+1)-old_x)*(A(n,2*k1+1)-new_x)*(new_y-old_y)...
                -(A(n-i,2*k1+1)-old_x)*(A(n-i,2*k1+1)-old_x)*(new_y-old_y)-(old_y-A(n-i,2*k1))*(A(n-i,2*k1+1)-old_x)*(new_x-old_x);
          end
          
        
        
        %d调整式子的分母
          d_denominator=(A(n,3)-new_x)^2+(A(n,5)-new_x)^2+(A(n,7)-new_x)^2+(A(n,9)-new_x)^2+(A(n,11)-new_x)^2+...
              (A(n-i,3)-old_x)^2+(A(n-i,5)-old_x)^2+(A(n-i,7)-old_x)^2+(A(n-i,9)-old_x)^2+(A(n-i,11)-old_x)^2;
       
        %利用d对最新的邻居消息进行调整（平移）
          d=d_molecule/d_denominator;
          for k2=1:5
              A(n,2*k2)=A(n,2*k2)+d;
          end
        end
     
        %线性回归求斜率
         Tj_mean=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10)+A(n-i,2)+A(n-i,4)+A(n-i,6)+A(n-i,8)+A(n-i,10))/10;
         Ti_mean=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11)+A(n-i,3)+A(n-i,5)+A(n-i,7)+A(n-i,9)+A(n-i,11))/10;
         
         molecule=0;
         for k3=1:5
             molecule=molecule+(A(n,2*k3)-Tj_mean)*(A(n,2*k3+1)-Ti_mean)+(A(n-i,2*k3)-Tj_mean)*(A(n-i,2*k3+1)-Ti_mean);
         end
         denominator=0;
         for k4=1:5
             denominator=denominator+(A(n,2*k4+1)-Ti_mean)^2+(A(n-i,2*k4+1)-Ti_mean)^2;
         end
    
       deltal(n-i,1)=molecule*A(n,12)/denominator;
%deltal(n-i,1)=(new_y-old_y+d_molecule/d_denominator)*A(n,12)/(new_x-old_x);
        m=m+1;

       %平均值窗口求斜率
%        alpha(1,1)=(A(n,10)-A(n-i,10))/(A(n,11)-A(n-i,11));
%        alpha(2,1)=(A(n,8)-A(n-i,8))/(A(n,9)-A(n-i,9));
%        alpha(3,1)=(A(n,6)-A(n-i,6))/(A(n,7)-A(n-i,7));
%        alpha(4,1)=(A(n,4)-A(n-i,4))/(A(n,5)-A(n-i,5));
%        alpha(5,1)=(A(n,2)-A(n-i,2))/(A(n,3)-A(n-i,3));
%        
%        deltal(n-i,1)=mean(alpha)*A(n,12);
%        m=m+1;

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
        deltal(n-i,1)=molecule*A(n-i,12)/denominator;
      
        m=m+1;
    end
end
l=(sum(deltal)+L)/(m+1);
end