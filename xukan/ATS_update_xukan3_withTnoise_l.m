function [l,ne]=ATS_update_xukan3_withTnoise_l(A,i,L,d)



    if(A(d+i,1)~=0)
        %对下n行消息进行偏移d的处理，让其尽量与第一次数据在同一直线上
        if A(d+i,1)-A(d,1)<=5 %当两次消息间隔轮次过大时，不需要进行偏移补偿
          old_x=(A(d,3)+A(d,5)+A(d,7)+A(d,9)+A(d,11))/5;%首次接收消息的x轴平均值
          old_y=(A(d,2)+A(d,4)+A(d,6)+A(d,8)+A(d,10))/5;%首次接收消息的y轴平均值
          new_x=(A(d+i,3)+A(d+i,5)+A(d+i,7)+A(d+i,9)+A(d+i,11))/5;%最新接收消息的x轴平均值
          new_y=(A(d+i,2)+A(d+i,4)+A(d+i,6)+A(d+i,8)+A(d+i,10))/5;%最新接收消息的y轴平均值
        
        %d调整的式子的分子
          d_molecule=0;
          for k1=1:5
              d_molecule=d_molecule+(A(d+i,2*k1)-old_y)*(new_x-old_x)*(A(d+i,2*k1+1)-...
                  new_x)-(A(d+i,2*k1+1)-old_x)*(A(d+i,2*k1+1)-new_x)*(new_y-old_y)...
                -(A(d,2*k1+1)-old_x)*(A(d,2*k1+1)-old_x)*(new_y-old_y)-...
                (old_y-A(d,2*k1))*(A(d,2*k1+1)-old_x)*(new_x-old_x);
          end
                  
        %d调整式子的分母
          d_denominator=(A(d+i,3)-new_x)^2+(A(d+i,5)-new_x)^2+(A(d+i,7)-new_x)^2+...
              (A(d+i,9)-new_x)^2+(A(d+i,11)-new_x)^2+...
              (A(d,3)-old_x)^2+(A(d,5)-old_x)^2+(A(d,7)-old_x)^2+...
              (A(d,9)-old_x)^2+(A(d,11)-old_x)^2;
       
        %利用d对最新的邻居消息进行调整（平移）
          d1=d_molecule/d_denominator;
          for k2=1:5
              A(d+i,2*k2)=A(d+i,2*k2)+d1;
          end
        end
     
        %线性回归求斜率
         Tj_mean=(A(d+i,2)+A(d+i,4)+A(d+i,6)+A(d+i,8)+A(d+i,10)+A(d,2)+A(d,4)+...
             A(d,6)+A(d,8)+A(d,10))/10;
         Ti_mean=(A(d+i,3)+A(d+i,5)+A(d+i,7)+A(d+i,9)+A(d+i,11)+A(d,3)+A(d,5)+...
             A(d,7)+A(d,9)+A(d,11))/10;
         
         molecule=0;
         for k3=1:5
             molecule=molecule+(A(d+i,2*k3)-Tj_mean)*(A(d+i,2*k3+1)-Ti_mean)+...
                 (A(d,2*k3)-Tj_mean)*(A(d,2*k3+1)-Ti_mean);
         end
         denominator=0;
         for k4=1:5
             denominator=denominator+(A(d+i,2*k4+1)-Ti_mean)^2+(A(d,2*k4+1)-Ti_mean)^2;
         end
         
       deltal=molecule*A(d+i,12)/denominator;%非首次接收的邻居的斜率，放入deltal中
       %将这个邻居信息加入返回值的邻居列表中
       ne=molecule*A(d+i,12)/denominator;
       if(A(d+i,14)~=0)
           l=(deltal+A(d+i,14)+L)/3;
       end
       if(A(d+i,14)==0)
           l=(deltal+L)/2;
       end
    end
    
    %第一次接收到该邻居消息且是最新消息
    if A(d+i,1)==0 && A(d,1)~=0
        %利用第一次的消息进行线性回归求斜率
       Tj_mean_old=(A(d,2)+A(d,4)+A(d,6)+A(d,8)+A(d,10))/5;
       Ti_mean_old=(A(d,3)+A(d,5)+A(d,7)+A(d,9)+A(d,11))/5;
        molecule=(A(d,2)-Tj_mean_old)*(A(d,3)-Ti_mean_old)+(A(d,4)-Tj_mean_old)*(A(d,5)-Ti_mean_old)...
            +(A(d,6)-Tj_mean_old)*(A(d,7)-Ti_mean_old)+(A(d,8)-Tj_mean_old)*(A(d,9)-Ti_mean_old)+...
            (A(d,10)-Tj_mean_old)*(A(d,11)-Ti_mean_old);
        denominator=(A(d,3)-Ti_mean_old)^2+(A(d,5)-Ti_mean_old)^2+(A(d,7)-Ti_mean_old)^2+...
            (A(d,9)-Ti_mean_old)^2+(A(d,11)-Ti_mean_old)^2;
        
        deltal=molecule*A(d,12)/denominator;%把首次接收到消息的节点的更新值加入deltal
        
        %将这个邻居信息加入返回值的邻居列表中
        ne=molecule*A(d,12)/denominator;
        if(A(d,14)~=0)
           l=(deltal+A(d,14)+L)/3;
       end
       if(A(d,14)==0)
           l=(deltal+L)/2;
       end
       
        
    end

end