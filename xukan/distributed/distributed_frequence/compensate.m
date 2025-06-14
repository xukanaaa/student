function [alphat,m1]=compensate(A,k,i,j)
%D1表示广播节点的用于同步的邻居信息矩阵
%k表示当前二跳节点是哪一个节点。
%i表示无人机总数
%j表示当前广播轮次
%返回值一个为斜率，另一个为计算斜率所用到的消息轮次

    
   %有第三轮该邻居最新消息
    if A(k+i+i,1)>=(j-1)
        %存在首轮消息，将二轮三轮与一轮距离补偿
        if A(k+i-i,1) ~=0
            
            %对第三轮消息进行偏移d的处理，让其尽量与第一次数据在同一直线上
            if abs(A(k+i+i,1)-A(k+i-i,1))<=10 %当两次消息间隔轮次过大时，不需要进行偏移补偿
                old_x=(A(k+i-i,3)+A(k+i-i,5)+A(k+i-i,7)+A(k+i-i,9)+A(k+i-i,11))/5;%首次接收消息的x轴平均值
                old_y=(A(k+i-i,2)+A(k+i-i,4)+A(k+i-i,6)+A(k+i-i,8)+A(k+i-i,10))/5;%首次接收消息的y轴平均值
                new_x=(A(k+i+i,3)+A(k+i+i,5)+A(k+i+i,7)+A(k+i+i,9)+A(k+i+i,11))/5;%最新接收消息的x轴平均值
                new_y=(A(k+i+i,2)+A(k+i+i,4)+A(k+i+i,6)+A(k+i+i,8)+A(k+i+i,10))/5;%最新接收消息的y轴平均值
                
                %d调整的式子的分子
                d_molecule=0;
                for k1=1:5
                    d_molecule=d_molecule+(A(k+i+i,2*k1)-old_y)*(new_x-old_x)*(A(k+i+i,2*k1+1)-...
                        new_x)-(A(k+i+i,2*k1+1)-old_x)*(A(k+i+i,2*k1+1)-new_x)*(new_y-old_y)...
                        -(A(k+i-i,2*k1+1)-old_x)*(A(k+i-i,2*k1+1)-old_x)*(new_y-old_y)-...
                        (old_y-A(k+i-i,2*k1))*(A(k+i-i,2*k1+1)-old_x)*(new_x-old_x);
                end
                
                %d调整式子的分母
                d_denominator=(A(k+i+i,3)-new_x)^2+(A(k+i+i,5)-new_x)^2+(A(k+i+i,7)-new_x)^2+...
                    (A(k+i+i,9)-new_x)^2+(A(k+i+i,11)-new_x)^2+...
                    (A(k+i-i,3)-old_x)^2+(A(k+i-i,5)-old_x)^2+(A(k+i-i,7)-old_x)^2+...
                    (A(k+i-i,9)-old_x)^2+(A(k+i-i,11)-old_x)^2;
                
                %利用d对最新的邻居消息进行调整（平移）
                d=d_molecule/d_denominator;
                for k2=1:5
                    A(k+i+i,2*k2)=A(k+i+i,2*k2)+d;
                end
            end
            
            %对第二轮消息进行偏移d的处理，让其尽量与第一次数据在同一直线上
            if abs(A(k+i,1)-A(k+i-i,1))<=10 %当两次消息间隔轮次过大时，不需要进行偏移补偿
                old_x=(A(k+i-i,3)+A(k+i-i,5)+A(k+i-i,7)+A(k+i-i,9)+A(k+i-i,11))/5;%首次接收消息的x轴平均值
                old_y=(A(k+i-i,2)+A(k+i-i,4)+A(k+i-i,6)+A(k+i-i,8)+A(k+i-i,10))/5;%首次接收消息的y轴平均值
                new_x=(A(k+i,3)+A(k+i,5)+A(k+i,7)+A(k+i,9)+A(k+i,11))/5;%最新接收消息的x轴平均值
                new_y=(A(k+i,2)+A(k+i,4)+A(k+i,6)+A(k+i,8)+A(k+i,10))/5;%最新接收消息的y轴平均值
                
                %d调整的式子的分子
                d_molecule=0;
                for k1=1:5
                    d_molecule=d_molecule+(A(k+i,2*k1)-old_y)*(new_x-old_x)*(A(k+i,2*k1+1)-...
                        new_x)-(A(k+i,2*k1+1)-old_x)*(A(k+i,2*k1+1)-new_x)*(new_y-old_y)...
                        -(A(k+i-i,2*k1+1)-old_x)*(A(k+i-i,2*k1+1)-old_x)*(new_y-old_y)-...
                        (old_y-A(k+i-i,2*k1))*(A(k+i-i,2*k1+1)-old_x)*(new_x-old_x);
                end
                
                %d调整式子的分母
                d_denominator=(A(k+i,3)-new_x)^2+(A(k+i,5)-new_x)^2+(A(k+i,7)-new_x)^2+...
                    (A(k+i,9)-new_x)^2+(A(k+i,11)-new_x)^2+...
                    (A(k+i-i,3)-old_x)^2+(A(k+i-i,5)-old_x)^2+(A(k+i-i,7)-old_x)^2+...
                    (A(k+i-i,9)-old_x)^2+(A(k+i-i,11)-old_x)^2;
                
                %利用d对最新的邻居消息进行调整（平移）
                d=d_molecule/d_denominator;
                for k2=1:5
                    A(k+i,2*k2)=A(k+i,2*k2)+d;
                end
            end
            
            %线性回归求斜率
            Tj_mean=(A(k+i,2)+A(k+i,4)+A(k+i,6)+A(k+i,8)+A(k+i,10)+A(k+i-i,2)+A(k+i-i,4)+...
                A(k+i-i,6)+A(k+i-i,8)+A(k+i-i,10)+A(k+i+i,2)+A(k+i+i,4)+A(k+i+i,6)+A(k+i+i,8)+A(k+i+i,10))/15;
            Ti_mean=(A(k+i,3)+A(k+i,5)+A(k+i,7)+A(k+i,9)+A(k+i,11)+A(k+i-i,3)+A(k+i-i,5)+...
                A(k+i-i,7)+A(k+i-i,9)+A(k+i-i,11)+A(k+i+i,3)+A(k+i+i,5)+A(k+i+i,7)+A(k+i+i,9)+A(k+i+i,11))/15;
            
            molecule=0;
            for k3=1:5
                molecule=molecule+(A(k+i,2*k3)-Tj_mean)*(A(k+i,2*k3+1)-Ti_mean)+...
                    (A(k+i-i,2*k3)-Tj_mean)*(A(k+i-i,2*k3+1)-Ti_mean)+(A(k+i+i,2*k3)-Tj_mean)*(A(k+i+i,2*k3+1)-Ti_mean);
            end
            denominator=0;
            for k4=1:5
                denominator=denominator+(A(k+i,2*k4+1)-Ti_mean)^2+(A(k+i-i,2*k4+1)-Ti_mean)^2+(A(k+i+i,2*k4+1)-Ti_mean)^2;
            end
            
            alphat=molecule/denominator;%非首次接收的邻居的斜率，放入deltal中
            m1=3;
            
            %首轮没有消息，表示后面两轮都是二跳消息，这时，不进行距离补偿
        else
            %线性回归求斜率
            Tj_mean=(A(k+i,2)+A(k+i,4)+A(k+i,6)+A(k+i,8)+A(k+i,10)+A(k+i+i,2)+A(k+i+i,4)+A(k+i+i,6)+A(k+i+i,8)+A(k+i+i,10))/10;
            Ti_mean=(A(k+i,3)+A(k+i,5)+A(k+i,7)+A(k+i,9)+A(k+i,11)+A(k+i+i,3)+A(k+i+i,5)+A(k+i+i,7)+A(k+i+i,9)+A(k+i+i,11))/10;
            
            molecule=0;
            for k3=1:5
                molecule=molecule+(A(k+i,2*k3)-Tj_mean)*(A(k+i,2*k3+1)-Ti_mean)+...
                    (A(k+i+i,2*k3)-Tj_mean)*(A(k+i+i,2*k3+1)-Ti_mean);
            end
            denominator=0;
            for k4=1:5
                denominator=denominator+(A(k+i,2*k4+1)-Ti_mean)^2+(A(k+i+i,2*k4+1)-Ti_mean)^2;
            end
            
            alphat=molecule/denominator;%非首次接收的邻居的斜率，放入deltal中
            m1=2;
        end
        
        
        %第二次接收到该邻居消息且是最新消息
    elseif  A(k+i,1)>=j-1
        %首轮有消息，对第一轮进行距离补偿
        if A(k+i-i,1)~=0
            %对第二轮消息进行偏移d的处理，让其尽量与第一次数据在同一直线上
            if abs(A(k+i,1)-A(k+i-i,1))<=10 %当两次消息间隔轮次过大时，不需要进行偏移补偿
                old_x=(A(k+i-i,3)+A(k+i-i,5)+A(k+i-i,7)+A(k+i-i,9)+A(k+i-i,11))/5;%首次接收消息的x轴平均值
                old_y=(A(k+i-i,2)+A(k+i-i,4)+A(k+i-i,6)+A(k+i-i,8)+A(k+i-i,10))/5;%首次接收消息的y轴平均值
                new_x=(A(k+i,3)+A(k+i,5)+A(k+i,7)+A(k+i,9)+A(k+i,11))/5;%最新接收消息的x轴平均值
                new_y=(A(k+i,2)+A(k+i,4)+A(k+i,6)+A(k+i,8)+A(k+i,10))/5;%最新接收消息的y轴平均值
                
                %d调整的式子的分子
                d_molecule=0;
                for k1=1:5
                    d_molecule=d_molecule+(A(k+i,2*k1)-old_y)*(new_x-old_x)*(A(k+i,2*k1+1)-...
                        new_x)-(A(k+i,2*k1+1)-old_x)*(A(k+i,2*k1+1)-new_x)*(new_y-old_y)...
                        -(A(k+i-i,2*k1+1)-old_x)*(A(k+i-i,2*k1+1)-old_x)*(new_y-old_y)-...
                        (old_y-A(k+i-i,2*k1))*(A(k+i-i,2*k1+1)-old_x)*(new_x-old_x);
                end
                
                %d调整式子的分母
                d_denominator=(A(k+i,3)-new_x)^2+(A(k+i,5)-new_x)^2+(A(k+i,7)-new_x)^2+...
                    (A(k+i,9)-new_x)^2+(A(k+i,11)-new_x)^2+...
                    (A(k+i-i,3)-old_x)^2+(A(k+i-i,5)-old_x)^2+(A(k+i-i,7)-old_x)^2+...
                    (A(k+i-i,9)-old_x)^2+(A(k+i-i,11)-old_x)^2;
                
                %利用d对最新的邻居消息进行调整（平移）
                d=d_molecule/d_denominator;
                for k2=1:5
                    A(k+i,2*k2)=A(k+i,2*k2)+d;
                end
            end
            
            %线性回归求斜率
            Tj_mean=(A(k+i,2)+A(k+i,4)+A(k+i,6)+A(k+i,8)+A(k+i,10)+A(k+i-i,2)+A(k+i-i,4)+...
                A(k+i-i,6)+A(k+i-i,8)+A(k+i-i,10))/10;
            Ti_mean=(A(k+i,3)+A(k+i,5)+A(k+i,7)+A(k+i,9)+A(k+i,11)+A(k+i-i,3)+A(k+i-i,5)+...
                A(k+i-i,7)+A(k+i-i,9)+A(k+i-i,11))/10;
            
            molecule=0;
            for k3=1:5
                molecule=molecule+(A(k+i,2*k3)-Tj_mean)*(A(k+i,2*k3+1)-Ti_mean)+...
                    (A(k+i-i,2*k3)-Tj_mean)*(A(k+i-i,2*k3+1)-Ti_mean);
            end
            denominator=0;
            for k4=1:5
                denominator=denominator+(A(k+i,2*k4+1)-Ti_mean)^2+(A(k+i-i,2*k4+1)-Ti_mean)^2;
            end
            
            alphat=molecule/denominator;%非首次接收的邻居的斜率，放入deltal中
            m1=2;
            %首轮没有消息，说明只有第二轮一个二跳信息，直接估计
        else
            %利用第一次的消息进行线性回归求斜率
        Tj_mean_old=(A(k+i,2)+A(k+i,4)+A(k+i,6)+A(k+i,8)+A(k+i,10))/5;
        Ti_mean_old=(A(k+i,3)+A(k+i,5)+A(k+i,7)+A(k+i,9)+A(k+i,11))/5;
        molecule=(A(k+i,2)-Tj_mean_old)*(A(k+i,3)-Ti_mean_old)+(A(k+i,4)-Tj_mean_old)*(A(k+i,5)-Ti_mean_old)...
            +(A(k+i,6)-Tj_mean_old)*(A(k+i,7)-Ti_mean_old)+(A(k+i,8)-Tj_mean_old)*(A(k+i,9)-Ti_mean_old)+...
            (A(k+i,10)-Tj_mean_old)*(A(k+i,11)-Ti_mean_old);
        denominator=(A(k+i,3)-Ti_mean_old)^2+(A(k+i,5)-Ti_mean_old)^2+(A(k+i,7)-Ti_mean_old)^2+...
            (A(k+i,9)-Ti_mean_old)^2+(A(k+i,11)-Ti_mean_old)^2;
        
        alphat=molecule*A(k+i,12)/denominator;%把首次接收到消息的节点的更新值加入deltal
        m1=1;
        end
        
        
        %第一次接收到该邻居消息且是最新消息
    elseif A(k+i-i,1)>=j-1
        %利用第一次的消息进行线性回归求斜率
        Tj_mean_old=(A(k+i-i,2)+A(k+i-i,4)+A(k+i-i,6)+A(k+i-i,8)+A(k+i-i,10))/5;
        Ti_mean_old=(A(k+i-i,3)+A(k+i-i,5)+A(k+i-i,7)+A(k+i-i,9)+A(k+i-i,11))/5;
        molecule=(A(k+i-i,2)-Tj_mean_old)*(A(k+i-i,3)-Ti_mean_old)+(A(k+i-i,4)-Tj_mean_old)*(A(k+i-i,5)-Ti_mean_old)...
            +(A(k+i-i,6)-Tj_mean_old)*(A(k+i-i,7)-Ti_mean_old)+(A(k+i-i,8)-Tj_mean_old)*(A(k+i-i,9)-Ti_mean_old)+...
            (A(k+i-i,10)-Tj_mean_old)*(A(k+i-i,11)-Ti_mean_old);
        denominator=(A(k+i-i,3)-Ti_mean_old)^2+(A(k+i-i,5)-Ti_mean_old)^2+(A(k+i-i,7)-Ti_mean_old)^2+...
            (A(k+i-i,9)-Ti_mean_old)^2+(A(k+i-i,11)-Ti_mean_old)^2;
        
        alphat=molecule*A(k+i-i,12)/denominator;%把首次接收到消息的节点的更新值加入deltal
        m1=1;
    end
end