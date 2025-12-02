%% 1. 初始化与参数
clear; clc; close all;

% === 通信参数 (已修改为 20MHz) ===
Fc = 20e6;                % 码片速率: 20 MHz (50 ns/chip)
Fs = 200e6;               % 采样率: 200 MHz (10倍过采样, 1点=5ns)
SpreadFactor = 255;       % Kasami 序列长度
SamplesPerChip = Fs/Fc;   % 10
SamplesPerSymbol = SpreadFactor * SamplesPerChip; 

% === 环境参数 ===
SNR_dB = 15;              % 信噪比
SIR_dB = -10;             % 强干扰
Fixed_Threshold = 0.5;    % 同步门限

% === 多径参数 (单位: 码片) ===
% Path 1: 主径
% Path 2: 0.75 码片 (37.5 ns) -> 造成主峰畸变
% Path 3: 3.0 码片 (150 ns)
% Path 4: 10.0 码片 (500 ns)
Multipath_Delays_Chips = [0,   0.75, 3.0,  10.0]; 
Multipath_Gains        = [1.0, 0.5,  0.3,  0.15]; 

% === 获取当前时间作为发送数据 ===
CurrentTime = datetime('now');
HH = CurrentTime.Hour;
MM = CurrentTime.Minute;
SS = round(CurrentTime.Second);
Tx_Time_Int = HH*10000 + MM*100 + SS; % 编码为整数

% 真实物理延迟
True_Delay_Samples = 1500.62; 

fprintf('=== 仿真开始 (20 MHz 码片速率) ===\n');
fprintf('发送时间数据: %02d:%02d:%02d\n', HH, MM, SS);
fprintf('真实延迟: %.2f Samples\n', True_Delay_Samples);

%% 2. 信号生成 (组帧)
% 1. 整数转二进制 (24位)
NumDataBits = 24; 
Tx_Bits = bitget(Tx_Time_Int, NumDataBits:-1:1)'; 
% 2. 加 Pilot
Tx_Frame = [1; Tx_Bits]; 
NumSymbols = length(Tx_Frame);

% 3. 扩频码生成
kasamiGen = comm.KasamiSequence('Polynomial', [8 4 3 2 0], 'SamplesPerFrame', SpreadFactor, 'InitialConditions', [0 0 0 0 0 0 0 1]);
kasamiGen.Index = 0; 
pnCode_Target = 2*kasamiGen() - 1; 

% 4. 调制
txSignal_Target = [];
for i = 1:NumSymbols
    if Tx_Frame(i) == 1, sym = pnCode_Target; else, sym = -pnCode_Target; end
    txSignal_Target = [txSignal_Target; rectpulse(sym, SamplesPerChip)]; %#ok<AGROW>
end
paddedTx = [zeros(500,1); txSignal_Target; zeros(500,1)];

% 5. 干扰生成
release(kasamiGen); kasamiGen.Index = 1; 
pnCode_Interf = 2*kasamiGen() - 1;
txSignal_Interf = [];
for i = 1:NumSymbols
    if rand > 0.5, sym = pnCode_Interf; else, sym = -pnCode_Interf; end
    txSignal_Interf = [txSignal_Interf; rectpulse(sym, SamplesPerChip)]; %#ok<AGROW>
end

refSignal = rectpulse(pnCode_Target, SamplesPerChip);

%% 3. 信道模拟 (多径 + 干扰)
sig_Target_Total = []; 
fDelayFilt = dsp.VariableFractionalDelay('InterpolationMethod', 'Farrow');

for k = 1:length(Multipath_Gains)
    Current_Delay = True_Delay_Samples + Multipath_Delays_Chips(k) * SamplesPerChip;
    frac = Current_Delay - floor(Current_Delay);
    int_d = floor(Current_Delay);
    
    sig_path = fDelayFilt(txSignal_Target, frac);
    sig_path = [zeros(int_d, 1); sig_path];
    
    if isempty(sig_Target_Total)
        sig_Target_Total = Multipath_Gains(k) * sig_path;
    else
        len_total = length(sig_Target_Total);
        len_path  = length(sig_path);
        if len_path > len_total
            sig_Target_Total = [sig_Target_Total; zeros(len_path - len_total, 1)];
        elseif len_total > len_path
            sig_path = [sig_path; zeros(len_total - len_path, 1)];
        end
        sig_Target_Total = sig_Target_Total + Multipath_Gains(k) * sig_path;
    end
end

Amp_Interf = sqrt(1 / 10^(SIR_dB/10)); 
sig_I = [zeros(500, 1); txSignal_Interf];
len = max(length(sig_Target_Total), length(sig_I)) + 500;
sig_Target_Total(end+1:len) = 0; sig_Target_Total = sig_Target_Total(1:len);
sig_I(end+1:len) = 0; sig_I = sig_I(1:len);

rxSignal = awgn(sig_Target_Total + Amp_Interf * sig_I, SNR_dB, 'measured');

%% 4. 接收处理
[acor_raw, lag] = xcorr(rxSignal, refSignal);
acor_norm = abs(acor_raw)/max(abs(acor_raw));

% 同步
[~, locs] = findpeaks(acor_norm, 'MinPeakHeight', Fixed_Threshold);
if isempty(locs), error('同步失败'); end
SyncIdx = locs(1); 

% TOA 估计 (Spline插值)
span = 6; 
idx_range = (SyncIdx-span):(SyncIdx+span);
lag_fine = linspace(lag(idx_range(1)), lag(idx_range(end)), 1000);
val_fine = interp1(lag(idx_range), acor_norm(idx_range), lag_fine, 'spline');
[~, fineMaxIdx] = max(val_fine);
Estimated_Delay = lag_fine(fineMaxIdx);
TOA_Error_ns = (Estimated_Delay - True_Delay_Samples)/Fs * 1e9;

% 解调
Rx_Bits = zeros(NumDataBits, 1);
for i = 1:NumDataBits
    sampIdx = SyncIdx + i * SamplesPerSymbol;
    if acor_raw(sampIdx) > 0, Rx_Bits(i) = 1; else, Rx_Bits(i) = 0; end
end
Rx_Time_Int = sum(Rx_Bits .* (2.^((NumDataBits-1):-1:0)'));

% 解析时间
Rx_HH = floor(Rx_Time_Int / 10000);
Rx_MM = floor(mod(Rx_Time_Int, 10000) / 100);
Rx_SS = mod(Rx_Time_Int, 100);

fprintf('--------------------------------------\n');
fprintf('解调时间: %02d:%02d:%02d\n', Rx_HH, Rx_MM, Rx_SS);
if Rx_Time_Int == Tx_Time_Int
    fprintf('状态: [成功]\n');
else
    fprintf('状态: [误码]\n');
end
fprintf('TOA 误差: %.4f ns\n', abs(TOA_Error_ns));
fprintf('--------------------------------------\n');

%% 5. 图形展示

% =======================================================
% Figure 1: 信号微观对比 (前 100 个码片)
% =======================================================
figure('Name', '1. Signal Micro View', 'Position', [100, 200, 800, 700]);
MicroChips = 100;
ZoomLen = MicroChips * SamplesPerChip;

subplot(2,1,1);
TxStart = 501; 
plot(1:ZoomLen, paddedTx(TxStart : TxStart + ZoomLen - 1), 'b', 'LineWidth', 1.2);
title(['(a) 发送数据包前段 (' sprintf('%02d:%02d:%02d', HH, MM, SS) ')']);
grid on; ylim([-1.5 1.5]); ylabel('幅度'); xlabel('局部采样点');
ax = gca; ax.XTick = 0 : 100 : ZoomLen; grid minor;

subplot(2,1,2);
RxStart = round(True_Delay_Samples); 
plot(1:ZoomLen, rxSignal(RxStart : RxStart + ZoomLen - 1), 'Color', [0.8 0.3 0.3]);
title(['(b) 接收信号微观视图 (20MHz, SIR=' num2str(SIR_dB) 'dB)']);
grid on; ylabel('幅度'); ylim([-6 6]); 
xlabel('局部采样点');
ax = gca; ax.XTick = 0 : 100 : ZoomLen; grid minor;

% =======================================================
% Figure 2: 互相关结果 (主峰原始形态 + 多径标记)
% =======================================================
figure('Name', '2. Correlation Peak (Focused)', 'Position', [950, 450, 800, 450]);

RangeSpan = 200; 
IdxZoom = (SyncIdx - RangeSpan) : (SyncIdx + RangeSpan);
CurrentLag = lag(IdxZoom);
CurrentVal = acor_raw(IdxZoom);

plot(CurrentLag, CurrentVal, 'b-', 'LineWidth', 1.5); hold on;
plot(lag(SyncIdx), acor_raw(SyncIdx), 'go', 'MarkerSize', 10, 'LineWidth', 2);
text(lag(SyncIdx), acor_raw(SyncIdx)+0.2, 'Sync Peak', 'Color', 'g', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% 标记多径位置 (仅虚线)
BaseLag = lag(SyncIdx);
MaxVal = max(CurrentVal);
xline(BaseLag + 0.75*SamplesPerChip, 'r--', 'LineWidth', 1.5); % Path 2
xline(BaseLag + 3.0*SamplesPerChip, 'm--', 'LineWidth', 1.5);  % Path 3
xline(BaseLag + 10.0*SamplesPerChip, 'm--', 'LineWidth', 1.5); % Path 4

yline(0, 'k-', 'Zero Line');
title(['互相关主峰 (时间 ' sprintf('%02d:%02d:%02d', Rx_HH, Rx_MM, Rx_SS) ' 解调基准)']);
ylabel('原始相关幅值'); xlabel('Lag (Samples)');
grid on; grid minor;
xlim([min(CurrentLag), max(CurrentLag)]);

% =======================================================
% Figure 3: 真实 vs 估计时间差距 (TOA Gap)
% =======================================================
figure('Name', '3. TOA Precision Gap', 'Position', [950, 50, 800, 350]);

idx_p = (SyncIdx-span):(SyncIdx+span);
lag_p = lag(idx_p);
val_p = acor_norm(idx_p);

t_axis = (lag_p - lag(SyncIdx))/Fs * 1e9;
t_fine = linspace(min(t_axis), max(t_axis), 500);
val_interp = interp1(t_axis, val_p, t_fine, 'spline');

plot(t_fine, val_interp, 'k-', 'LineWidth', 1.5); hold on;
stem(t_axis, val_p, 'b-o', 'LineWidth', 1.2, 'MarkerFaceColor', 'b');

t_true_rel = (True_Delay_Samples - lag(SyncIdx))/Fs * 1e9;
t_est_rel  = (Estimated_Delay - lag(SyncIdx))/Fs * 1e9;

xline(t_true_rel, 'g-', 'LineWidth', 2);
xline(t_est_rel, 'r--', 'LineWidth', 2);

y_mid = max(val_p) * 0.96;
quiver(t_true_rel, y_mid, t_est_rel-t_true_rel, 0, 0, 'Color', 'm', 'LineWidth', 2, 'MaxHeadSize', 0.5);
text((t_true_rel+t_est_rel)/2, y_mid+0.005, sprintf('Gap: %.3f ns', abs(TOA_Error_ns)), ...
    'Color', 'm', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

title('TOA 估计精度 (20MHz, 纳秒级)');
xlabel('相对时间 (ns)'); ylabel('归一化幅度');
grid on; xlim([-3, 3]); ylim([0.9, 1.01]);