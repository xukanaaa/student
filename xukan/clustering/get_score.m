%% 计算节点作为簇头的得分
function score=get_score(i,A,num,j)
% i:自身节点编号
% A:自身的连接矩阵
% num:无人机总数
% j:当前轮次
% score：节点得分

%计算自己的频率稳定度
score_fre=0;
if i>=1&&i<=5
    score_fre=5;
elseif i>=6&&i<=15
    score_fre=10;
elseif i>=16&&i<=30
    score_fre=15;
else
    score_fre=20;
end

%计算节点的度
dn=0;
for ii=1:num
    if A(i,ii)==1
        dn=dn+1;
    end
end

%计算加权得分
score=score_fre-dn/100;
end