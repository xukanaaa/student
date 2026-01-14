function [x_k, P_k, is_jump] = clock_sync_kf(x_km1, P_km1, dt, t_send, t_recv)
% CLOCK_SYNC_KF_ADAPTIVE 自适应卡尔曼滤波时间同步算法
%
% 改进点：去除了 'round' 参数，改为通过残差检测自动识别频偏跳变
%
% 输入:
%   x_km1  : 上一时刻状态 [Offset (s); Skew (s/s)]
%   P_km1  : 上一时刻协方差
%   dt     : 时间间隔 (s)
%   t_send : 发送时间
%   t_recv : 接收时间
%
% 输出:
%   x_k    : 更新状态
%   P_k    : 更新协方差
%   is_jump: 调试标志，1表示本轮检测到了突变

    % --- 1. 基础参数设置 ---
    % 测量噪声 R (5ns 的平方)
    R = 2.5e-17; 
    
    % --- 2. 定义两套过程噪声参数 ---
    % 方案 A: 平稳模式 (Steady Mode)
    % 物理含义：高稳晶振，几乎无漂移。
    Sf_steady = 1e-23;
    Sg_steady = 1e-26; % 极小，保证直线非常平滑
    
    % 方案 B: 突变模式 (Jump Mode)
    % 物理含义：用于应对 0.2ppm (2e-7) 级别的瞬间跳变
    % 这里的 Sg 需要足够大，让 P_skew 瞬间膨胀到 1e-12 ~ 1e-13 级别
    Sf_jump = 1e-16;
    Sg_jump = 1e-12;   % 关键参数：允许频率瞬间发生 1ppm 级的变化
    
    % 构建平稳状态下的 Q 矩阵
    Q_steady = [Sf_steady*dt + Sg_steady*(dt^3)/3, Sg_steady*(dt^2)/2; 
                Sg_steady*(dt^2)/2,                Sg_steady*dt];

    % 构建突变状态下的 Q 增量矩阵 (或者直接用它来膨胀 P)
    Q_jump_add = [Sf_jump*dt + Sg_jump*(dt^3)/3, Sg_jump*(dt^2)/2; 
                  Sg_jump*(dt^2)/2,              Sg_jump*dt];

    % --- 3. 状态转移矩阵 F ---
    F = [1, dt; 
         0, 1];
    
    % --- 4. 预测阶段 (Predict) ---
    x_pred = F * x_km1;
    
    % 默认先按“平稳模式”进行预测
    P_pred = F * P_km1 * F' + Q_steady;
    
    % --- 5. 观测与残差检测 (Detection) ---
    H = [1, 0];
    z = t_recv - t_send;
    
    % 初步计算残差 (Innovation)
    y = z - H * x_pred;
    
    % 初步计算理论残差方差 S
    S_cov = H * P_pred * H' + R;
    
    % [关键改进] 计算归一化残差平方 (NIS)
    % NIS 服从卡方分布。对于标量测量，NIS = y^2 / S
    % 如果模型准确，NIS 应该很小（<9.0, 即 3-sigma 范围）
    nis = (y^2) / S_cov;
    
    % 设定检测阈值
    % 3-sigma -> 9.0;  4-sigma -> 16.0;  5-sigma -> 25.0
    % 建议设为 16.0 或 25.0 以避免噪声误触发
    nis_threshold = 20.0; 
    
    is_jump = 0;
    
    % --- 6. 自适应调整 (Adaptation) ---
    if nis > nis_threshold
        % 触发条件：残差太大，说明发生了突变，原先的 P_pred 太小了（太自信了）
        is_jump = 1;
        
        % 策略：强行“膨胀”预测协方差矩阵 P_pred
        % 这告诉滤波器：“我现在很不确定状态，请更多地相信当前的测量值 z”
        P_pred = P_pred + Q_jump_add;
        
        % [重要] 因为 P_pred 变了，必须重新计算 S_cov 和卡尔曼增益 K
        S_cov = H * P_pred * H' + R;
        
        % (可选) 重新计算残差？不需要，因为 x_pred 没变，y 没变
    end
    
    % --- 7. 更新状态 (Update) ---
    % 计算卡尔曼增益 K
    K = P_pred * H' / S_cov;
    
    % 更新状态向量
    x_k = x_pred + K * y;
    
    % 更新协方差矩阵
    I = eye(2);
    P_k = (I - K * H) * P_pred;
    
    % --- 8. 强制对称 ---
    P_k = (P_k + P_k') / 2;
    
end



% function [x_k, P_k] = clock_sync_kf(x_km1, P_km1, dt, t_send, t_recv)
% % CLOCK_SYNC_KF 卡尔曼滤波时间同步算法 (秒级原单位版)
% %
% % 输入:
% %   x_km1  : 上一时刻状态 [2x1] -> [Offset (s); Skew (s/s)]
% %            x(1) 是时间偏差，x(2) 是频率偏差(如 1e-6 代表 1ppm)
% %   P_km1  : 上一时刻协方差 [2x2]
% %   dt     : 距离上一时刻的时间间隔 (s)
% %   t_send : 发送时间 (s)
% %   t_recv : 接收时间 (s) (已扣除传播延迟)
% %
% % 输出:
% %   x_k    : 更新后的状态 [Offset (s); Skew (s/s)]
% %   P_k    : 更新后的协方差
% 
%     % --- 1. 参数设置 (标准单位: 秒) ---
%     % 过程噪声强度 (Process Noise)
%     % 物理含义：晶振频率的不稳定性。
%     % 普通晶振约为 1e-9 到 1e-8 量级，平方后为 1e-16 到 1e-18
%     S = 0; 
% 
%     % 测量噪声方差 (Measurement Noise R)
%     % 物理含义：时间戳打标误差的平方。
%     % 假设误差为 10ns (1e-8 s)，则 R = 1e-16
%     R = 2.5e-17; 
% 
%     % --- 2. 状态转移矩阵 F ---
%     % Offset_k = Offset_{k-1} + Skew_{k-1} * dt
%     % Skew_k   = Skew_{k-1}
%     F = [1, dt; 
%          0, 1];
% 
%     % --- 3. 过程噪声矩阵 Q ---
%     % 采用白噪声加速模型 (White Noise Acceleration)
%     Q = [(1e-23)*dt+(1e-26)*((dt)^3)/3, (1e-26)*((dt)^2)/2; 
%              (1e-26)*((dt)^2)/2, (1e-26)*dt];
% 
%     % --- 4. 预测阶段 (Predict) ---
%     x_pred = F * x_km1;
%     P_pred = F * P_km1 * F' + Q;
% 
%     % --- 5. 观测阶段 (Update) ---
%     % 观测矩阵 H: 我们观测到的是时间偏差 offset
%     H = [1, 0];
% 
%     % 观测值 z: 当前测量到的时间偏差 (接收 - 发送)
%     z = t_recv - t_send;
% 
%     % 计算残差 (Innovation)
%     y = z - H * x_pred;
% 
%     % 计算残差协方差 S_cov
%     S_cov = H * P_pred * H' + R;
% 
%     % --- 6. 数值稳定性保护 ---
%     % 因为单位是秒，S_cov 可能非常小 (e.g. 1e-15)。
%     % 但如果它小于 MATLAB 的 eps 或接近 0，除法会溢出。
%     % if abs(S_cov) < 1e-20
%     %     % 测量无效或方差过小，保持预测值
%     %     x_k = x_pred;
%     %     P_k = P_pred;
%     %     return;
%     % end
% 
%     % --- 7. 更新状态 ---
%     % 卡尔曼增益 K
%     K = P_pred * H' / S_cov;
% 
%     % 更新状态向量
%     x_k = x_pred + K * y;
% 
%     % 更新协方差矩阵 (Joseph form 保证正定性，或简单形式)
%     I = eye(2);
%     P_k = (I - K * H) * P_pred;
% 
%     % --- 8. 强制对称 ---
%     % 防止多次迭代后浮点误差导致矩阵不对称
%     P_k = (P_k + P_k') / 2;
% 
% end