function h=gtsp_update_h(A,j,i,H)
%GTSP的逻辑时钟相位调整
delta_h=zeros(i,1);
m=1;
for n=1:i
    if A(n,1)>=(j-1)&&A(n,6)~=0
        delta_h(n,1)=A(n,4)-A(n,5);
        m=m+1;
    end
end
h=H+(sum(delta_h))/m;
end