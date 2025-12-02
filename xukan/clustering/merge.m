%% 合并自身和邻居的分数矩阵，赋值给邻居
function S=merge(S_neb,S_self,num)
% S_neb:邻居的分数矩阵
% S_self:自己的分数矩阵
% num:无人机总数
% S:合并后的矩阵
S=zeros(num,1);
for i=1:num
    if S_neb(i,1)>0
        S(i,1)=S_neb(i,1);
    end
    if S_self(i,1)>0
        S(i,1)=S_self(i,1);
    end
end

end