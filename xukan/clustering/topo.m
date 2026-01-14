function draw_timing_diagram_fixed()
    % 创建画布
    figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    hold on;
    axis off; 
    
    % --- 参数设置 ---
    y_node_i = 6;  % Node i 的高度
    y_node_j = 2;  % Node j 的高度
    x_limit = [0, 12]; % 时间轴长度
    
    % --- 1. 画两条时间轴 ---
    plot(x_limit, [y_node_i, y_node_i], 'k-', 'LineWidth', 1.5);
    plot(x_limit, [y_node_j, y_node_j], 'k-', 'LineWidth', 1.5);
    
    % 轴标签
    text(0.5, y_node_i + 0.5, 'Node $i$', 'Interpreter', 'latex', 'FontSize', 14);
    text(0.5, y_node_j - 0.5, 'Node $j$', 'Interpreter', 'latex', 'FontSize', 14);

    % --- 2. 画左侧的垂直关系 ---
    x_left = 2; 
    draw_arrow([x_left, x_left], [y_node_i, y_node_j], 'k'); 
    text(x_left - 0.2, (y_node_i + y_node_j)/2, 'level $i$ + $T_{send}$', ...
        'Interpreter', 'latex', 'HorizontalAlignment', 'right');
    
    plot([x_left+0.5, x_left+1.5], [y_node_j, y_node_i], 'k-');
    draw_arrow_head(x_left+1.5, y_node_i, x_left+0.5, y_node_j, 0.3);
    text(x_left+0.8, (y_node_i + y_node_j)/2, 'level $j$', 'Interpreter', 'latex');

    % --- 3. 画中间的 Delta t 区域 ---
    t_ref = 4.5; 
    t_shift = 5.0; 
    
    plot([t_ref, t_ref], [y_node_i, y_node_j-1], 'k--'); 
    plot([t_shift, t_shift], [y_node_i-1, y_node_j], 'k--'); 
    
    y_dt = y_node_j - 0.8;
    draw_double_arrow([t_ref, t_shift], [y_dt, y_dt], 0.1);
    text((t_ref+t_shift)/2, y_dt - 0.4, '$\Delta t$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');

    % --- 4. 画交叉传输 ---
    Ti_1 = 5.0;  
    Ti_2 = 9.0;  
    Tj_1 = 7.0;  
    Tj_send = 5.8; 
    
    % Node i -> Node j
    plot([Ti_1, Tj_1], [y_node_i, y_node_j], 'k-', 'LineWidth', 1.2);
    draw_arrow_head(Tj_1, y_node_j, Ti_1, y_node_i, 0.4); 
    text(Ti_1, y_node_i + 0.3, '$T_i(1)$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
    text(Tj_1, y_node_j - 0.3, '$T_j(1)$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
    
    % Node j -> Node i
    plot([Tj_send, Ti_2], [y_node_j, y_node_i], 'k-', 'LineWidth', 1.2);
    draw_arrow_head(Ti_2, y_node_i, Tj_send, y_node_j, 0.4); 
    text(Ti_2, y_node_i + 0.3, '$T_i(2)$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
    text(Tj_send + 1.5, y_node_j - 0.3, '$T_j(2)$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');

    % --- 5. 画大括号 (已修正函数) ---
    % 上方大括号
    draw_curly_brace(Ti_1, Ti_2, y_node_i + 0.6, 0.3, 1); 
    text((Ti_1 + Ti_2)/2, y_node_i + 1.2, '$\rho_j^i$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % 下方大括号
    draw_curly_brace(Tj_send, Tj_1, y_node_j - 0.6, 0.3, -1); 
    text((Tj_send + Tj_1)/2, y_node_j - 1.2, '$\rho_i^j$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 12);

    ylim([0, 8]);
    xlim([0, 11]);
    set(gca, 'Visible', 'off'); 
end

% --- 辅助函数 ---

function draw_arrow(x, y, color)
    quiver(x(1), y(1), x(2)-x(1), y(2)-y(1), 0, 'Color', color, 'MaxHeadSize', 0.1, 'LineWidth', 1.2);
end

function draw_arrow_head(x_tip, y_tip, x_tail, y_tail, scale)
    theta = atan2(y_tip - y_tail, x_tip - x_tail);
    alpha = pi/8; 
    r = scale; 
    x1 = x_tip - r * cos(theta - alpha);
    y1 = y_tip - r * sin(theta - alpha);
    x2 = x_tip - r * cos(theta + alpha);
    y2 = y_tip - r * sin(theta + alpha);
    patch([x_tip, x1, x2], [y_tip, y1, y2], 'k');
end

function draw_double_arrow(x, y, h)
    plot(x, y, 'k-');
    draw_arrow_head(x(1), y(1), x(2), y(2), h*3);
    draw_arrow_head(x(2), y(2), x(1), y(1), h*3);
end

% --- 修正后的画大括号函数 ---
function draw_curly_brace(x1, x2, y, h, dir)
    % x1, x2: 起点和终点 X
    % y: 基准高度
    % h: 大括号的高度幅度
    % dir: 1 (向上), -1 (向下)
    
    mid = (x1 + x2) / 2;
    w = x2 - x1;
    
    % 修正点：确保 xx 中的值是唯一的且递增
    % 我们不使用重复点，而是用很小的偏移量 (w*0.05) 来制造"转角"
    xx = [x1, x1+w*0.05, mid-w*0.05, mid, mid+w*0.05, x2-w*0.05, x2];
    
    % 对应的 Y 坐标形状：低-高-高-尖峰-高-高-低
    % h_base 是括号肩膀的高度，h_tip 是中间尖尖的高度
    h_base = h * 0.7; 
    h_tip  = h * 1.0; 
    
    yy = [y, y+dir*h_base, y+dir*h_base, y+dir*h_tip, y+dir*h_base, y+dir*h_base, y];
    
    % 生成平滑曲线
    xx_smooth = linspace(x1, x2, 100);
    % 使用 'pchip' 插值，因为它比 spline 更适合保持形状（不会过冲）
    yy_smooth = interp1(xx, yy, xx_smooth, 'pchip');
    
    plot(xx_smooth, yy_smooth, 'k-', 'LineWidth', 1);
end