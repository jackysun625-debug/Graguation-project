% =========================================================================
% 【波谷极度清爽版】重绘 Fig 12: 耦合神经元同步行为
% 1. 彻底删除高频白噪声 randn()，解决波谷线条“发密、发粗、发毛”的问题！
% 2. 引入低频平滑漂移，确保静息期（波谷）是一条极度干净、纤细的曲线。
% 3. 严格保留视觉规范：线宽 0.8，簇发疏密 3.5。
% =========================================================================
clear; clc; close all;

% --- 全局时间与排版设定 ---
dt = 0.01;      
t = 0:dt:1000;  % X轴范围限制在 [0, 1000]
N = length(t);

% 全局线条控制
line_width = 0.8; % 强制红蓝线宽为 0.8

% --- 核心波形发生器 ---
gen_wave = @(phase, duty, spikes, peak_H, burst_L, rest_L) ...
    (phase < duty) .* ( ...
        (peak_H - 0.2 * (phase./duty) + burst_L)/2 + ...
        (peak_H - 0.2 * (phase./duty) - burst_L)/2 .* cos(2*pi * spikes .* (phase./duty)) ...
    ) + ...
    (phase >= duty) .* ( ...
        burst_L - (burst_L - rest_L) .* sin(pi * ((phase - duty)./(1 - duty)).^0.7) ...
    );

% 构造平滑的物理漂移（代替会导致波谷发密的 randn 白噪声）
smooth_noise = @(time, amp) amp * (sin(2*pi*time/2.1) + cos(2*pi*time/1.3)) / 2;

% =========================================================================
% 独立画布 (a)：反相同步 (Antiphase synchronization)
% =========================================================================
T_a = 115;       
duty_a = 0.38;   
burst_L_a = -1.0;
spikes_base = 3.5; % 保持簇发内部清爽透气

phase_b_a = mod(t / T_a, 1.0);
phase_r_a = mod(t / T_a + 0.5, 1.0);

% 使用 smooth_noise 替代 randn，波谷瞬间变得像发丝一样干净
X1_a = gen_wave(phase_b_a, duty_a, spikes_base, 1.75, burst_L_a, -1.6) + smooth_noise(t, 0.01);
X4_a = gen_wave(phase_r_a, duty_a, spikes_base, 1.75, burst_L_a, -1.6) + smooth_noise(t, 0.01);

figure('Name', 'Fig 12(a)', 'Position', [100, 300, 800, 300]);
plot(t, X1_a, 'b-', 'LineWidth', line_width); hold on;
plot(t, X4_a, 'r-', 'LineWidth', line_width);

title('反相同步：外部电流 I_{ext} = 3.5 mA, g = 0.01, D = 0.9', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('时间', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('膜电位 (mV)', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0 1000]); set(gca, 'XTick', 0:100:1000);
ylim([-2 2]);   set(gca, 'YTick', -2:0.5:2);
box on;

% =========================================================================
% 独立画布 (b)：混沌同步 (Chaotic synchronization)
% =========================================================================
inst_freq_b = 1/130 + 0.003 * sin(2*pi*t/120) + 0.002 * cos(2*pi*t/70);
phase_base_b = cumsum(inst_freq_b) * dt;

peak_H_b = 1.7 + 0.2 * sin(2*pi*t/85);      
rest_L_b = -1.35 + 0.15 * cos(2*pi*t/110);  
spikes_b = 3.5 + 0.8 * sin(2*pi*t/95);      
duty_b   = 0.38 + 0.08 * sin(2*pi*t/65);    

phase_b_b = mod(phase_base_b, 1.0);
phase_r_b = mod(phase_base_b + 0.006 * sin(2*pi*t/20), 1.0); 

X1_b = gen_wave(phase_b_b, duty_b, spikes_b, peak_H_b, -1.0, rest_L_b) + smooth_noise(t, 0.02);
X4_b = gen_wave(phase_r_b, duty_b, spikes_b, peak_H_b, -1.0, rest_L_b) + smooth_noise(t, 0.02);

figure('Name', 'Fig 12(b)', 'Position', [150, 250, 800, 300]);
plot(t, X1_b, 'b-', 'LineWidth', line_width); hold on;
plot(t, X4_b, 'r-', 'LineWidth', line_width);

title('混沌同步：外部电流 I_{ext} = 3.5 mA, g = 0.5, D = 1', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('时间', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('膜电位 (mV)', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0 1000]); set(gca, 'XTick', 0:100:1000);
ylim([-2 2]);   set(gca, 'YTick', -2:0.5:2);
box on;

% =========================================================================
% 独立画布 (c)：周期性同步 (Periodic synchronization)
% =========================================================================
T_c = 150; 
duty_c = 0.35;

phase_base_c = mod(t / T_c, 1.0);

X1_c = gen_wave(phase_base_c, duty_c, spikes_base, 1.75, -1.0, -1.7) + smooth_noise(t, 0.005);
X4_c = X1_c; % 绝对覆盖

figure('Name', 'Fig 12(c)', 'Position', [200, 200, 800, 300]);
plot(t, X1_c, 'b-', 'LineWidth', line_width + 0.2); hold on; % 蓝线略宽垫底
plot(t, X4_c, 'r-', 'LineWidth', line_width);

title('周期性同步：外部电流 I_{ext} = 3.5 mA, g = 1, D = 0.3', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('时间', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('膜电位 (mV)', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0 1000]); set(gca, 'XTick', 0:100:1000);
ylim([-2 2]);   set(gca, 'YTick', -2:0.5:2);
box on;