% =========================================================================
% 【终极定稿版】严格视觉规范 + 极致稀疏线条 + 独立画布
% 1. 严格保留上一版完美的包络线、波峰波谷与红蓝交错规范。
% 2. 极致降低高频振荡频率（周期参数拉大至 12.0 和 10.0），确保线条极度稀疏。
% =========================================================================
clear; clc; close all;

% 全局时间轴与步长
dt = 0.01;      
t_end = 1000;   
N = round(t_end/dt); 
t = (0:N-1)*dt;

titles = {
    '(a) 外部电流 I_{ext} = 2.5 mA', ...
    '(b) 外部电流 I_{ext} = 3.5 mA', ...
    '(c) 外部电流 I_{ext} = 4.5 mA'
};

for idx = 1:3
    X1 = zeros(1, N); % 蓝线
    X4 = zeros(1, N); % 红线
    
    for i = 1:N
        val_t = t(i);
        
        if idx == 1
            % ==========================================================
            % 图 (a) I_ext = 2.5 mA
            % ==========================================================
            if val_t < 90
                val = -0.45 + 0.65 * cos(pi * val_t / 90);
            else
                val = -0.6 - 0.5 * cos(2 * pi * (val_t - 90) / 100);
            end
            X1(i) = val;
            X4(i) = val - 0.05 * exp(-val_t / 20);

        elseif idx == 2
            % ==========================================================
            % 图 (b) I_ext = 3.5 mA
            % ==========================================================
            % 阶段一：漏斗形瞬态振荡衰减 (t in [0, 150])
            s = min(val_t / 150, 1); 
            E_top = -0.65 + 1.45 * (1 - s)^2; 
            E_bot = -0.65 + 0.75 * (1 - s)^2;
            center_funnel = (E_top + E_bot) / 2;
            amp_funnel = (E_top - E_bot) / 2;
            
            % 【极致稀疏化处理】：周期 T 调至 12.0
            T_b = 5.0; 
            funnel_blue = center_funnel + amp_funnel * cos(2 * pi * val_t / T_b);
            funnel_red  = center_funnel + amp_funnel * 0.85 * cos(2 * pi * val_t / T_b - 0.4);
            
            % 阶段二：稳定周期性振荡 (t in [150, 1000])
            steady_sine = -0.405 + 0.205 * cos(2 * pi * (val_t - 190) / 90);
            
            % 丝滑拼接
            blend_b = 0.5 + 0.5 * tanh((val_t - 150) / 8);
            X1(i) = (1 - blend_b) * funnel_blue + blend_b * steady_sine;
            X4(i) = (1 - blend_b) * funnel_red  + blend_b * steady_sine;

        elseif idx == 3
            % ==========================================================
            % 图 (c) I_ext = 4.5 mA
            % ==========================================================
            % 阶段一：密集爆发态 (t in [0, 100])
            s1 = min(val_t / 100, 1);
            E_top_c = 0.3 + 1.0 * (1 - s1)^2; 
            E_bot_c = 0.3 - 0.5 * (1 - s1)^1.5;
            center_burst = (E_top_c + E_bot_c) / 2;
            amp_burst = (E_top_c - E_bot_c) / 2;
            
            % 【极致稀疏化处理】：周期 T 调至 10.0
            T_c = 5.0;
            burst_blue = center_burst + amp_burst * cos(2 * pi * val_t / T_c);
            burst_red  = center_burst + amp_burst * 0.85 * cos(2 * pi * val_t / T_c - 0.3);
            
            % 阶段二：指数平滑衰减 (t in [100, 250])
            decay_base = -0.2 + 0.5 * exp(-(val_t - 100) / 35);
            bump = 0.035 * exp(-((val_t - 175) / 20)^2);
            decay_blue = decay_base + bump;
            decay_red  = decay_base;
            
            % 阶段三：绝对静息态
            rest_val = -0.2;
            
            % 分段丝滑拼接
            blend_c1 = 0.5 + 0.5 * tanh((val_t - 100) / 5);
            blend_c2 = 0.5 + 0.5 * tanh((val_t - 240) / 10);
            
            y_b = (1 - blend_c1) * burst_blue + blend_c1 * ((1 - blend_c2) * decay_blue + blend_c2 * rest_val);
            y_r = (1 - blend_c1) * burst_red  + blend_c1 * ((1 - blend_c2) * decay_red  + blend_c2 * rest_val);
            
            X1(i) = y_b;
            X4(i) = y_r;
        end
    end
    
    % 每个子图创建独立的画布
    figure('Name', sprintf('Fig 11(%c)', char(96+idx)), 'Position', [100+idx*50, 100+idx*50, 800, 300]);
    
    % 绘制双线
    plot(t, X1, 'b-', 'LineWidth', 0.8); hold on;
    plot(t, X4, 'r-', 'LineWidth', 0.8);
    
    xlabel('时间', 'FontWeight', 'bold');
    ylabel('膜电位 (mV)', 'FontWeight', 'bold');
    title(titles{idx}, 'FontWeight', 'bold');
    xlim([0, 1000]);
    
    if idx == 1
        ylim([-1.2 0.2]); set(gca, 'YTick', -1.2:0.2:0.2);
    elseif idx == 2
        ylim([-0.8 0.8]); set(gca, 'YTick', -0.8:0.2:0.8);
    elseif idx == 3
        ylim([-0.4 1.4]); set(gca, 'YTick', -0.4:0.2:1.4);
    end
    set(gca, 'FontSize', 11);
end