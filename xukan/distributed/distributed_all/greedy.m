%在当前状态以贪婪策略选择下次的动作
function action_curr_RL=greedy(state_curr_RL,Qtable)
% state_curr_RL：当前的状态索引
% Qtable：Q表格

vals = [ -5, -2, -1, 0, 1, 2, 5];
%生成随机数，百分之95概率大于一
p = rand;
if p < 0.95
    % 90%概率：生成大于1的随机数（例如1到10之间）
    epsilon = 1 + 9*rand;  % 范围(1,10)
else
    % 10%概率：生成小于1的随机数（例如0到1之间）
    epsilon = rand;        % 范围(0,1)
end

if epsilon>1
    % 提取目标行
    row_data = Qtable(state_curr_RL, :);

    % 找到最大值
    max_val = max(row_data);

    % 找到所有最大值的索引
    max_indices = find(row_data == max_val);

    % 随机选择一个索引（若只有一个则直接返回）
    random_index = max_indices(randi(length(max_indices)));
    action_curr_RL=vals(random_index);
else
    action_curr_RL=vals(randi([1, 7], 1));
end
end