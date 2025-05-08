function l=MulBro_OneRouWithoutTnoise_l(A,j,i,L)
%GTSP的逻辑时钟频率调整
deltal=zeros(i,1);%每个节点对应的频率调整值
alpha=zeros(4,1);%未滤波时的每个节点对应的十个相对频率偏移
m=0;
for n=1:i
    %相邻点作差求平均值计算斜率
    if (A(n,1)>=(j-1))
        for k1=1:4
            alpha(k1)=(A(n,2*k1+2)-A(n,2*k1))/(A(n,2*k1+3)-A(n,2*k1+1))+2e-6*randn;
        end
        deltal(n,1)=mean(alpha)*A(n,12);
        m=m+1;
    end
end
l=(sum(deltal)+L)/(m+1);
end