function M=MTS_l_h(Fre,Pha,j,num,L,H,i)
% Fre:频率消息矩阵
% Pha:相位消息矩阵
% j:当前同步轮次
% num:无人机总数
% L:自身节点频率补偿值
% H:自身节点相位补偿值
% i:自身节点编号
% M:返回值，(1,1)表示l补偿值，(2,1)表示h补偿值,(1,3)表示同步状态

%初始化返回值矩阵
M=zeros(1,3);
L_old=L;
H_old=H;


%遍历矩阵Fre,进行频率偏移补偿值计算和相位偏移补偿值计算
for n=num+1:2*num
    %有多于两轮同一节点消息，才用于同步
    if Fre(n,1)>=j-1
        thre=(9e-5)/(Fre(n,1)-Fre(n-num,1));

        alpha_ij=(Fre(n,2)-Fre(n-num,2))/(Fre(n,3)-Fre(n-num,3));
        q_ij=alpha_ij*Pha(n-num,4)/L;
        %邻居的补偿后的频偏大于自己的，更新自己的频偏与时钟值为邻居的
        if q_ij>1+thre
            M(1,1)=alpha_ij*Pha(n-num,4);
            M(1,2)=Pha(n-num,4)*Pha(n-num,2)+Pha(n-num,5)-alpha_ij*Pha(n-num,4)*Pha(n-num,3);
        end
        %二者频偏一致，将自己时钟调整为二者最大值
        if q_ij>=1-thre&&q_ij<=1+thre
            M(1,1)=L;
            M(1,2)=max(L*Pha(n-num,3)+H,Pha(n-num,4)*Pha(n-num,2)+Pha(n-num,5))-L*Pha(n-num,3);
        end
        %自己的频偏比邻居大，自身的不变
        if q_ij<1-thre
            M(1,1)=L;
            M(1,2)=H;
        end
        L=M(1,1);
        H=M(1,2);
    end
end

if L_old==L&&abs(H_old-H)<6e-6
    if Fre(i,5)==0
        M(1,3)=0.5;
    end
    if Fre(i,5)==0.5
        M(1,3)=1;
    end
    if Fre(i,5)==1
        M(1,3)=1;
    end
else
    M(1,3)=0;
end

end