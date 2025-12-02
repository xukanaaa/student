% 定义X轴数据（3000到7500，每500一个点）
x = 3000:500:7500;

% 定义Y轴数据
data1 = [21, 31, 41, 56, 66, 86, 111, 131, 151, 176];  % 数据1
data2 = [21, 26, 31, 36, 41, 46, 56, 66, 81, 121];   % 数据2

% 创建图形
figure;
hold on; grid on;  % 保持图形并添加网格

% 绘制两条曲线
plot(x, data1, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(x, data2, 's-', 'LineWidth', 1.5, 'MarkerSize', 6);

% 设置坐标轴范围和标签
xlim([3000, 7500]);
xlabel('网络直径（米）', 'FontSize', 10);
ylabel('同步收敛需要的同步轮数', 'FontSize', 10);

% 设置标题
title('网络同步收敛速度仿真', 'FontSize', 12);

% 添加图例
legend('数据1', '数据2', 'Location', 'best');

% 优化显示
box on;  % 给坐标轴添加边框
hold off;