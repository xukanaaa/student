%从相位矩阵获取方差
function v=get_var(A,num,L,H)
% A：强化学习阶段用于同步的相位矩阵
% num：无人机总数
M=zeros(num,1);%初始化时钟差值矩阵
B=zeros(2,5);%第一行表示邻居的逻辑时钟，第二行表示自身的逻辑时钟
%遍历矩阵A
for n=1:num
    if A(n,1)==1
        
            A(n,2)=A(n,2)-A(n,14);%传播时延补偿
        
        
            B(1,1)=A(n,2)*A(n,12)+A(n,13);
            B(2,1)=A(n,3)*L+H;
        
        M(n,1)=B(1,1)-B(2,1);
    end
end

if all(all(M == 0))
    v=0;
else
v=var(M(M ~= 0));
end

end