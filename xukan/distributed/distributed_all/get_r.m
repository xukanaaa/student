%从动作，方差，转移到的状态三个维度计算价值r
function r=get_r(action,var_curr_RL,state_curr_RL)
% action:上一个状态采取了什么动作
% var_curr_RL：采取这个动作后下次同步时的方差
% state_curr_RL：下一个状态
v_max=5e-10;%允许的最大方差
%精度维度参数
r_k1=10;
r_alpha=1e6;
r_beta=1e7;
r_c=50;
%周期维度参数
r_gamma=0.9;
r_v0=1e-10;
r_epsilon=1;
%趋势维度参数
r_D=5;
r_E=1e5;
if var_curr_RL<=v_max
    %精度维度的奖励
    r_accuracy=r_k1-r_alpha*var_curr_RL;
    %周期维度奖励
    r_T=r_gamma*action*(r_v0/(var_curr_RL+r_epsilon));
    %趋势维度奖励
    if mod(state_curr_RL,3)==1
        r_trend=-1*r_E*var_curr_RL;
    elseif mod(state_curr_RL,3)==2
        r_trend=r_D;
    else
        r_trend=r_D/2;
    end
else
    r_accuracy=-1*r_beta*(var_curr_RL-v_max)-r_c;
    r_T=0;
    %趋势维度奖励
    if mod(state_curr_RL,3)==1
        r_trend=-1*r_E*var_curr_RL;
    elseif mod(state_curr_RL,3)==2
        r_trend=r_D;
    else
        r_trend=r_D/2;
    end
end
r=r_accuracy+r_T+r_trend;
end