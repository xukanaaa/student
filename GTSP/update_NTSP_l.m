function l=update1_l(A,j,i,L)
%基于邻居选择的逻辑时钟频率调整
if j==2
   deltal=zeros(i,1);
   m=0;
   for n=1:i
       if (A(n,1)>=(j-1))&&(A(n,6)~=0)
           deltal(n,1)=((A(n,2)-A(n,6))*A(n,8))/(A(n,3)-A(n,7));
           m=m+1;
       end
   end
   l=(sum(deltal)+L)/(m+1);
end

if j>2
    syn=0;
    unsyn=0;
    delta_syn=zeros(i,1);
    delta_unsyn=zeros(i,1);
    for k=1:i
        if (A(k,1)>=(j-1))&&(A(k,6)~=0)
            delta_syn(k,1)=((A(k,2)-A(k,6))*A(k,8))/(A(k,3)-A(k,7));
            syn=syn+1;
        end
        if (A(k,1)<=(-j+1))&&(A(k,6)~=0)
            delta_unsyn(k,1)=((A(k,2)-A(k,6))*A(k,8))/(A(k,3)-A(k,7));
            unsyn=unsyn+1;
        end
    end
    if syn>=(2*unsyn)
        l=(sum(delta_syn)+L)/(syn+1);
    end
    if syn<(2*unsyn)
        l=(sum(delta_syn)+sum(delta_unsyn)+L)/(syn+unsyn+1);
    end
end