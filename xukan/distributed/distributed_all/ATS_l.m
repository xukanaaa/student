function l=ATS_l(A,j,L)
%ATS的逻辑时钟频率调整
%A：收到广播节点的邻居信息矩阵
%j：当前广播节点的编号
%L：收到广播节点自身的逻辑时钟频率调整参数
deltal=L;
    if (A(j,6)~=0)
        deltal=((A(j,2)-A(j,6))/(A(j,3)-A(j,7)))*A(j,8);
    end
l=(deltal+L)/2;
end
