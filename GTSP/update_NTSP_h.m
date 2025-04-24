function h=update1_h(A,j,i,H)
%基于邻居选择的逻辑时钟相位调整
if j==2
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

if j>2
    syn=0;
    unsyn=0;
    delta_syn=zeros(i,1);
    delta_unsyn=zeros(i,1);
    for k=1:i
        if (A(k,1)>=(j-1))&&(A(k,6)~=0)
            delta_syn(k,1)=A(k,4)-A(k,5);
            syn=syn+1;
        end
        if (A(k,1)<=(-j+1))&&(A(k,6)~=0)
            delta_unsyn(k,1)=A(k,4)-A(k,5);
            unsyn=unsyn+1;
        end
    end
    if syn>=(2*unsyn)
       h=H+(sum(delta_syn))/(syn+1);
    end
    if syn<(2*unsyn)
        h=H+(sum(delta_syn)+sum(delta_unsyn))/(syn+unsyn+1);
    end
end