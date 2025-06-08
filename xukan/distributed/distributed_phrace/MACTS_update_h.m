function h=MACTS_update_h(A,j,i,H,comm_range)
%GTSP的逻辑时钟相位调整
deltah=zeros(i,1);
m=1;
for n=1:i
    if A(n,1)>=(j-1)&&A(n,6)~=0
        deltah(n,1)=A(n,4)-A(n,5)-(comm_range/2)/3e8;
        m=m+1;
    end
end
h=H+(sum(deltah))/m;
end