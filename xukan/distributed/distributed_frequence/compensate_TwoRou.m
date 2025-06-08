function [alphat,m1]=compensate_TwoRou(A,k,i,j)
%A表示广播节点的用于同步的邻居信息矩阵
%k表示当前二跳节点是哪一个节点。
%i表示无人机总数
%j表示当前广播轮次
%返回值一个为斜率，另一个为计算斜率所用到的消息轮次

n=k+i;
%有第二轮该邻居最新消息（二跳，一跳）（二跳，空）（一跳，一跳）
if A(n,1)>=(j-1)
    %一跳，一跳
    if A(n,14)==0
        %对第二轮消息进行偏移d的处理，让其尽量与第一次数据在同一直线上
        if A(n,1)-A(n-i,1)<20 %当两次消息间隔轮次过大时，不需要进行偏移补偿
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
        
        alphat=molecule/denominator;%非首次接收的邻居的斜率，放入deltal中
        m1=2;
    end
    %最新的是虚拟，第一轮没有消息
    if A(n,14)>0&&A(n-i,1)==0
        %线性回归求斜率
        Tj_mean=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10))/5;
        Ti_mean=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11))/5;
        
        molecule=0;
        for k3=1:5
            molecule=molecule+(A(n,2*k3)-Tj_mean)*(A(n,2*k3+1)-Ti_mean);
        end
        denominator=0;
        for k4=1:5
            denominator=denominator+(A(n,2*k4+1)-Ti_mean)^2;
        end
        
        alphat=molecule/denominator;%非首次接收的邻居的斜率，放入deltal中
        m1=1;
    end
    %第二轮是虚拟，且存在第一轮消息(不是最新)，不进行距离补偿，使用两轮消息同步
    if A(n,14)>0&&A(n-i,1)>0&&A(n-i,1)<j-1
        if A(n,1)-A(n-i,1)<20 %当两次消息间隔轮次过大时，不需要进行偏移补偿
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
        
        alphat=molecule/denominator;%非首次接收的邻居的斜率，放入deltal中
        m1=2;
    end
    %第二轮是虚拟，且第一轮是最新消息，使用第一轮消息同步
    if A(n,14)>0&&A(n-i,1)>=j-1
        %利用第一次的消息进行线性回归求斜率
        Tj_mean_old=(A(n-i,2)+A(n-i,4)+A(n-i,6)+A(n-i,8)+A(n-i,10))/5;
        Ti_mean_old=(A(n-i,3)+A(n-i,5)+A(n-i,7)+A(n-i,9)+A(n-i,11))/5;
        molecule=(A(n-i,2)-Tj_mean_old)*(A(n-i,3)-Ti_mean_old)+(A(n-i,4)-Tj_mean_old)*(A(n-i,5)-Ti_mean_old)...
            +(A(n-i,6)-Tj_mean_old)*(A(n-i,7)-Ti_mean_old)+(A(n-i,8)-Tj_mean_old)*(A(n-i,9)-Ti_mean_old)+...
            (A(n-i,10)-Tj_mean_old)*(A(n-i,11)-Ti_mean_old);
        denominator=(A(n-i,3)-Ti_mean_old)^2+(A(n-i,5)-Ti_mean_old)^2+(A(n-i,7)-Ti_mean_old)^2+...
            (A(n-i,9)-Ti_mean_old)^2+(A(n-i,11)-Ti_mean_old)^2;
        
        alphat=molecule/denominator;%把首次接收到消息的节点的更新值加入deltal
        m1=1;
    end
    
else
    %第二轮没有最新消息（真实或者虚拟），且第一轮有最新消息
    if A(n-i,1)>=(j-1)
        %利用第一次的消息进行线性回归求斜率
        Tj_mean_old=(A(n-i,2)+A(n-i,4)+A(n-i,6)+A(n-i,8)+A(n-i,10))/5;
        Ti_mean_old=(A(n-i,3)+A(n-i,5)+A(n-i,7)+A(n-i,9)+A(n-i,11))/5;
        molecule=(A(n-i,2)-Tj_mean_old)*(A(n-i,3)-Ti_mean_old)+(A(n-i,4)-Tj_mean_old)*(A(n-i,5)-Ti_mean_old)...
            +(A(n-i,6)-Tj_mean_old)*(A(n-i,7)-Ti_mean_old)+(A(n-i,8)-Tj_mean_old)*(A(n-i,9)-Ti_mean_old)+...
            (A(n-i,10)-Tj_mean_old)*(A(n-i,11)-Ti_mean_old);
        denominator=(A(n-i,3)-Ti_mean_old)^2+(A(n-i,5)-Ti_mean_old)^2+(A(n-i,7)-Ti_mean_old)^2+...
            (A(n-i,9)-Ti_mean_old)^2+(A(n-i,11)-Ti_mean_old)^2;
        
        alphat=molecule/denominator;%把首次接收到消息的节点的更新值加入deltal
        m1=1;
    end
end
end