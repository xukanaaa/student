function h=update_h(A,j,i,H)
%GTSP的逻辑时钟相位调整
deltah=zeros(i,1);
m=1;
for n=1:i
    if (A(n,1)>=(j-1)&&(A(n,6)~=0))
        deltah(n,1)=A(n,4)-A(n,5);
        m=m+1;
    end
end
h=H+(sum(deltah))/m;
end