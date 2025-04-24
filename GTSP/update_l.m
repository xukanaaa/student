function l=update_l(A,j,i,L)
%GTSP的逻辑时钟频率调整
deltal=zeros(i,1);
m=0;
for n=1:i
    if (A(n,1)>=(j-1))&&(A(n,6)~=0)
        deltal(n,1)=((A(n,2)-A(n,6))/(A(n,3)-A(n,7))+2e-6*randn)*A(n,8);
        m=m+1;
    end
end
l=(sum(deltal)+L)/(m+1);
end
