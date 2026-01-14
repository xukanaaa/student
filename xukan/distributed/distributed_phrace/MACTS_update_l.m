function l=macts_update_l(A,j,i,L)
%GTSP的逻辑时钟频率调整
delta_l=zeros(i,1);
m=0;
for n=1:i
    if (A(n,1)>=(j-1))&&(A(n,6)~=0)
        delta_l(n,1)=((A(n,2)-A(n,6))/(A(n,3)-A(n,7)))*A(n,8);
        m=m+1;
    end
end
l=(sum(delta_l)+L)/(m+1);
end
