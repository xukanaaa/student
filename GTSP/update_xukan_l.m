function l=update_xukan_l(A,j,i,L)
%GTSP的逻辑时钟频率调整
deltal=zeros(i,1);%每个节点对应的频率调整值
alpha=zeros(4,1);%未滤波时的每个节点对应的十个相对频率偏移
m=0;
for n=1:i
    if (A(n,1)>=(j-1))
        %Tj_mean=(A(n,2)+A(n,4)+A(n,6)+A(n,8)+A(n,10))/5;
        %Ti_mean=(A(n,3)+A(n,5)+A(n,7)+A(n,9)+A(n,11))/5;
        %molecule=(A(n,2)-Tj_mean)*(A(n,3)-Ti_mean)+(A(n,4)-Tj_mean)*(A(n,5)-Ti_mean)...
         %   +(A(n,6)-Tj_mean)*(A(n,7)-Ti_mean)+(A(n,8)-Tj_mean)*(A(n,9)-Ti_mean)+...
          %  (A(n,10)-Tj_mean)*(A(n,11)-Ti_mean);
        %denominator=(A(n,3)-Ti_mean)^2+(A(n,5)-Ti_mean)^2+(A(n,7)-Ti_mean)^2+...
         %   (A(n,9)-Ti_mean)^2+(A(n,11)-Ti_mean)^2;
        %deltal(n,1)=molecule*A(n,12)/denominator;
        alpha(1)=(A(n,4)-A(n,2))/(A(n,5)-A(n,3))+2e-6*randn;
        %alpha(2)=(A(n,6)-A(n,2))/(A(n,7)-A(n,3));%+2e-6*randn;
        %alpha(3)=(A(n,8)-A(n,2))/(A(n,9)-A(n,3));%+2e-6*randn;
        %alpha(4)=(A(n,10)-A(n,2))/(A(n,11)-A(n,3));%+2e-6*randn;
        alpha(2)=(A(n,6)-A(n,4))/(A(n,7)-A(n,5))+2e-6*randn;
        %alpha(6)=(A(n,8)-A(n,4))/(A(n,9)-A(n,5));%+2e-6*randn;
        %alpha(7)=(A(n,10)-A(n,4))/(A(n,11)-A(n,5));%+2e-6*randn;
        alpha(3)=(A(n,8)-A(n,6))/(A(n,9)-A(n,7))+2e-6*randn;
        %alpha(9)=(A(n,10)-A(n,6))/(A(n,11)-A(n,7));%+2e-6*randn;
        alpha(4)=(A(n,10)-A(n,8))/(A(n,11)-A(n,9))+2e-6*randn;
        deltal(n,1)=mean(alpha)*A(n,12);
        m=m+1;
    end
end
l=(sum(deltal)+L)/(m+1);
end