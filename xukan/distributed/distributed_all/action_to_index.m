%从动作转化为索引
function action_index=action_to_index(action)
% action:动作，比如-5，-2，-1,0,1,2,5
map = containers.Map([-5,-2,-1,0,1,2,5], 1:7);
action_index = map(action); 
end