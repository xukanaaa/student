function [alphat,m1]=compensate(A,k,i,j)
%D1表示广播节点的用于同步的邻居信息矩阵
%k表示当前二跳节点是哪一个节点。
%i表示无人机总数
%j表示当前广播轮次
%返回值一个为斜率，另一个为计算斜率所用到的消息轮次

    
    %有第三轮该邻居最新消息
    if A(2*i+k,1)>=(j-1)
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
        Tj_mean=(A(k+i,2)+A(k+i,4)+A(k+i,6)+A(k+i,8)+A(k+i,10)+A(k,2)+A(k,4)+...
            A(k,6)+A(k,8)+A(k,10)+A(k+i*2,2)+A(k+i*2,4)+A(k+i*2,6)+A(k+i*2,8)+A(k+i*2,10))/15;
        Ti_mean=(A(k+i,3)+A(k+i,5)+A(k+i,7)+A(k+i,9)+A(k+i,11)+A(k,3)+A(k,5)+...
            A(k,7)+A(k,9)+A(k,11)+A(k+i*2,3)+A(k+i*2,5)+A(k+i*2,7)+A(k+i*2,9)+A(k+i*2,11))/15;
        
        molecule=0;
        for k3=1:5
            molecule=molecule+(A(k+i,2*k3)-Tj_mean)*(A(k+i,2*k3+1)-Ti_mean)+...
                (A(k,2*k3)-Tj_mean)*(A(k,2*k3+1)-Ti_mean)+(A(k+i*2,2*k3)-Tj_mean)*(A(k+i*2,2*k3+1)-Ti_mean);
        end
        denominator=0;
        for k4=1:5
            denominator=denominator+(A(k+i,2*k4+1)-Ti_mean)^2+(A(k,2*k4+1)-Ti_mean)^2+(A(k+i*2,2*k4+1)-Ti_mean)^2;
        end
        
        alphat=molecule/denominator;%非首次接收的邻居的斜率，放入deltal中
        m1=3;
    
    
    %第二次接收到该邻居消息且是最新消息
    elseif  A(i+k,1)>=j-1
        %线性回归求斜率
        Tj_mean=(A(k+i,2)+A(k+i,4)+A(k+i,6)+A(k+i,8)+A(k+i,10)+A(k,2)+A(k,4)+...
            A(k,6)+A(k,8)+A(k,10))/10;
        Ti_mean=(A(k+i,3)+A(k+i,5)+A(k+i,7)+A(k+i,9)+A(k+i,11)+A(k,3)+A(k,5)+...
            A(k,7)+A(k,9)+A(k,11))/10;
        
        molecule=0;
        for k3=1:5
            molecule=molecule+(A(k+i,2*k3)-Tj_mean)*(A(k+i,2*k3+1)-Ti_mean)+...
                (A(k,2*k3)-Tj_mean)*(A(k,2*k3+1)-Ti_mean);
        end
        denominator=0;
        for k4=1:5
            denominator=denominator+(A(k+i,2*k4+1)-Ti_mean)^2+(A(k,2*k4+1)-Ti_mean)^2;
        end
        
        alphat=molecule/denominator;%非首次接收的邻居的斜率，放入deltal中
        m1=2;
    
    
    %第一次接收到该邻居消息且是最新消息
    elseif A(k,1)>=j-1
             %利用第一次的消息进行线性回归求斜率
            Tj_mean_old=(A(k,2)+A(k,4)+A(k,6)+A(k,8)+A(k,10))/5;
            Ti_mean_old=(A(k,3)+A(k,5)+A(k,7)+A(k,9)+A(k,11))/5;
             molecule=(A(k,2)-Tj_mean_old)*(A(k,3)-Ti_mean_old)+(A(k,4)-Tj_mean_old)*(A(k,5)-Ti_mean_old)...
                 +(A(k,6)-Tj_mean_old)*(A(k,7)-Ti_mean_old)+(A(k,8)-Tj_mean_old)*(A(k,9)-Ti_mean_old)+...
                 (A(k,10)-Tj_mean_old)*(A(k,11)-Ti_mean_old);
             denominator=(A(k,3)-Ti_mean_old)^2+(A(k,5)-Ti_mean_old)^2+(A(k,7)-Ti_mean_old)^2+...
                 (A(k,9)-Ti_mean_old)^2+(A(k,11)-Ti_mean_old)^2;
    
             alphat=molecule/denominator;%把首次接收到消息的节点的更新值加入deltal
             m1=1;
    end
end