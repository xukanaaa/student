function l=update_xukan1_withTnoise_l(A,j,i,L)
%GTSP的逻辑时钟频率调整
deltal=zeros(i,1);%每个节点对应的频率调整值
%alpha=zeros(4,1);%未滤波时的每个节点对应的十个相对频率偏移
m=0;
for n=1:i
    %线性回归求散点斜率
    if (A(n,1)>=(j-1))
        Tj_mean=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10))/5;
        Ti_mean=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11))/5;
        molecule=0;
        for k1=1:5
            molecule=molecule+(A(n,2*k1)-Tj_mean)*(A(n,2*k1+1)-Ti_mean);
        end
        denominator=0;
        for k1=1:5
            denominator=denominator+(A(n,2*k1+1)-Ti_mean)^2;
        end
        deltal(n,1)=molecule*A(n,12)/denominator;
        m=m+1;
    end
end
l=(sum(deltal)+L)/(m+1);
end