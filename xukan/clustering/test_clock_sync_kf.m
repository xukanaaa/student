function [x_k, P_k, nis] = test_clock_sync_kf(x_km1, P_km1, dt, t_send, t_recv)
% ROBUST_CLOCK_SYNC_KF 抗差自适应时钟同步算法
%
% 针对场景：
% 1. 高稳晶振：日常运行极平稳
% 2. 温漂跟踪：允许 1ppm 级别的频偏变化 (Adaptive)
% 3. 鲁棒性：抵抗异常数据注入或网络风暴 (Robustness)
%
% 输入:
%   x_km1      : 上一时刻状态 [Offset (s); Skew (s/s)]
%   P_km1      : 上一时刻协方差
%   dt         : 时间间隔
%   t_send, t_recv : 时间戳
%
% 输出:
%   x_k, P_k   : 更新后的状态与协方差
%   status_code: 0=正常, 1=温漂/突变(已适应), 2=异常/攻击(已拒绝)

    % --- 1. 参数定义 ---
    % 测量噪声标准差 5ns -> 方差 R
    sigma_meas = 5e-9;
    R = sigma_meas^2; 
    
    % 过程噪声 (Q) 参数
    % 平稳模式 (高稳晶振)
    Sf_steady = 1e-23; 
    Sg_steady = 1e-26; % 极低，让 Skew 几乎不变
    
    % 适应模式 (用于跟踪 1ppm 温漂)
    % 1ppm = 1e-6。若 dt=1s，Skew 变化 1e-6 导致 Offset 误差累积。
    % 我们需要释放对 Skew 的信任，允许它快速变化。
    Sf_adapt = 1e-16;
    Sg_adapt = 1e-12;  % 允许频偏快速调整
    
    % 构建基础 Q 矩阵 (稳态)
    Q_steady = [Sf_steady*dt + Sg_steady*(dt^3)/3, Sg_steady*(dt^2)/2; 
                Sg_steady*(dt^2)/2,                Sg_steady*dt];
                
    % 构建适应性 Q 增量 (用于膨胀)
    Q_adapt_add = [Sf_adapt*dt + Sg_adapt*(dt^3)/3, Sg_adapt*(dt^2)/2; 
                   Sg_adapt*(dt^2)/2,              Sg_adapt*dt];

    % --- 2. 预测 (Predict) ---
    F = [1, dt; 
         0, 1];
         
    x_pred = F * x_km1;
    P_pred = F * P_km1 * F' + Q_steady; % 默认按稳态预测
    
    % --- 3. 观测与 NIS 计算 ---
    H = [1, 0];
    z = t_recv - t_send;
    
    y = z - H * x_pred;          % 新息 (Innovation)
    S = H * P_pred * H' + R;     % 新息协方差
    nis = (y^2) / S;             % 归一化新息平方 (NIS)
    
    % --- 4. 双阈值判决逻辑 (核心修改) ---
    % 阈值设定基于卡方分布 (Chi-squared, 1 degree of freedom)
    % Threshold_Low:  约 4-sigma (NIS=16.0)，超过此值认为发生温漂/跳变
    % Threshold_High: 约 10-sigma (NIS=100.0)，超过此值认为是非法攻击/野值
    
    thresh_adapt = 20;  
    thresh_reject = 5e8; % 根据实际网络抖动情况，这个值可能需要调大到 400-900
    
    status_code = 0; % 默认正常
    
    % 逻辑分支
    if nis > thresh_reject
        % [Case A: 拒绝] 异常数据或攻击
        % 策略：完全不更新状态，保持预测值（或者仅做极微小的更新）
        status_code = 2;
        
        x_k = x_pred;
        P_k = P_pred;
        
        % 直接返回，防止野值污染系统
        return; 
        
    elseif nis > thresh_adapt
        % [Case B: 适应] 温漂或真实跳变
        % 策略：膨胀协方差，强迫滤波器信任当前测量值
        status_code = 1;
        
        % 膨胀 P_pred，重点是增加对 Skew 变化的不确定性
        P_pred = P_pred + Q_adapt_add;
        
        % 重新计算 S (因为 P_pred 变了)
        S = H * P_pred * H' + R;
        
        % 继续向下执行更新...
    end
    
    % --- 5. 状态更新 (Update) ---
    % 计算卡尔曼增益
    K = P_pred * H' / S;
    
    % 更新状态
    x_k = x_pred + K * y;
    
    % 更新协方差 - 使用 Joseph Form (数值稳定性更高)
    % P = (I-KH)P(I-KH)' + KRK'
    % 在 R 极小 (5e-9) 的情况下，标准公式容易导致 P 非正定
    I = eye(2);
    ImKH = I - K * H;
    P_k = ImKH * P_pred * ImKH' + K * R * K';
    
    % 强制对称
    P_k = (P_k + P_k') / 2;
    
end