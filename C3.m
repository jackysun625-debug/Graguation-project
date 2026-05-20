% =========================================================================
% 改进H-R模型 ISI（峰间距）分岔图 
% 严格匹配文献 Fig. 3: Bifurcation diagram of ISI versus current
% =========================================================================
clear; clc; close all;

% 1. 模型全局参数设置 (依据原文献)
a = 1; b = 3; c = 1; d = 5; s = 4; r = 0.006; x0 = -1.6;
k = 0.4; k1 = 1; k2 = 0.5; alpha = 0.1; beta = 0.02; gamma = 0.2;

% 2. 仿真时间与步长设置
dt = 0.02;          % 积分步长
t_transient = 1500; % 舍弃前1500s瞬态，确保系统进入稳态
t_observe = 2500;   % 观察期长度，要足够长以捕捉大ISI
N_trans = round(t_transient/dt);
N_obs = round(t_observe/dt);

% 3. 扫描参数设置
I_range = 1.0 : 0.02 : 4.5; % 扫描范围匹配文献Fig. 3的X轴
bif_I = [];                 % 用于保存散点图的X坐标(外部电流)
bif_ISI = [];               % 用于保存散点图的Y坐标(ISI峰间距)

% 创建进度条
h = waitbar(0, '正在极速计算分岔图，请稍候...');
fprintf('开始计算分岔图，预计需要几十秒...\n');

% 4. 开始遍历外部电流 (这就是报错提示的第 28 行)
for idx = 1:length(I_range)
    I_ext = I_range(idx);
    
    % 初始化状态变量 [x1, x2, x3, phi]
    x1 = -1.5; x2 = 0; x3 = 0; phi = 0;
    
    % 第一阶段：跑掉瞬态过程 (将RK4公式直接内联展开，大幅提高MATLAB运算速度)
    for i = 1:N_trans
        % k1
        rho1 = alpha*phi^2 + beta*phi + gamma;
        dx1_1 = x2 - a*x1^3 + b*x1^2 - x3 + I_ext - k*rho1*x1;
        dx2_1 = c - d*x1^2 - x2;
        dx3_1 = r*(s*(x1 - x0) - x3);
        dphi_1 = k1*x1 - k2*phi;
        
        % k2
        x1_2 = x1 + 0.5*dt*dx1_1; x2_2 = x2 + 0.5*dt*dx2_1; x3_2 = x3 + 0.5*dt*dx3_1; phi_2 = phi + 0.5*dt*dphi_1;
        rho2 = alpha*phi_2^2 + beta*phi_2 + gamma;
        dx1_2 = x2_2 - a*x1_2^3 + b*x1_2^2 - x3_2 + I_ext - k*rho2*x1_2;
        dx2_2 = c - d*x1_2^2 - x2_2;
        dx3_2 = r*(s*(x1_2 - x0) - x3_2);
        dphi_2 = k1*x1_2 - k2*phi_2;
        
        % k3
        x1_3 = x1 + 0.5*dt*dx1_2; x2_3 = x2 + 0.5*dt*dx2_2; x3_3 = x3 + 0.5*dt*dx3_2; phi_3 = phi + 0.5*dt*dphi_2;
        rho3 = alpha*phi_3^2 + beta*phi_3 + gamma;
        dx1_3 = x2_3 - a*x1_3^3 + b*x1_3^2 - x3_3 + I_ext - k*rho3*x1_3;
        dx2_3 = c - d*x1_3^2 - x2_3;
        dx3_3 = r*(s*(x1_3 - x0) - x3_3);
        dphi_3 = k1*x1_3 - k2*phi_3;
        
        % k4
        x1_4 = x1 + dt*dx1_3; x2_4 = x2 + dt*dx2_3; x3_4 = x3 + dt*dx3_3; phi_4 = phi + dt*dphi_3;
        rho4 = alpha*phi_4^2 + beta*phi_4 + gamma;
        dx1_4 = x2_4 - a*x1_4^3 + b*x1_4^2 - x3_4 + I_ext - k*rho4*x1_4;
        dx2_4 = c - d*x1_4^2 - x2_4;
        dx3_4 = r*(s*(x1_4 - x0) - x3_4);
        dphi_4 = k1*x1_4 - k2*phi_4;
        
        % 更新状态
        x1 = x1 + (dt/6)*(dx1_1 + 2*dx1_2 + 2*dx1_3 + dx1_4);
        x2 = x2 + (dt/6)*(dx2_1 + 2*dx2_2 + 2*dx2_3 + dx2_4);
        x3 = x3 + (dt/6)*(dx3_1 + 2*dx3_2 + 2*dx3_3 + dx3_4);
        phi = phi + (dt/6)*(dphi_1 + 2*dphi_2 + 2*dphi_3 + dphi_4);
    end
    
    % 第二阶段：记录稳态观察期数据
    x1_obs = zeros(1, N_obs);
    for i = 1:N_obs
        rho1 = alpha*phi^2 + beta*phi + gamma; dx1_1 = x2 - a*x1^3 + b*x1^2 - x3 + I_ext - k*rho1*x1; dx2_1 = c - d*x1^2 - x2; dx3_1 = r*(s*(x1 - x0) - x3); dphi_1 = k1*x1 - k2*phi;
        x1_2 = x1 + 0.5*dt*dx1_1; x2_2 = x2 + 0.5*dt*dx2_1; x3_2 = x3 + 0.5*dt*dx3_1; phi_2 = phi + 0.5*dt*dphi_1; rho2 = alpha*phi_2^2 + beta*phi_2 + gamma; dx1_2 = x2_2 - a*x1_2^3 + b*x1_2^2 - x3_2 + I_ext - k*rho2*x1_2; dx2_2 = c - d*x1_2^2 - x2_2; dx3_2 = r*(s*(x1_2 - x0) - x3_2); dphi_2 = k1*x1_2 - k2*phi_2;
        x1_3 = x1 + 0.5*dt*dx1_2; x2_3 = x2 + 0.5*dt*dx2_2; x3_3 = x3 + 0.5*dt*dx3_2; phi_3 = phi + 0.5*dt*dphi_2; rho3 = alpha*phi_3^2 + beta*phi_3 + gamma; dx1_3 = x2_3 - a*x1_3^3 + b*x1_3^2 - x3_3 + I_ext - k*rho3*x1_3; dx2_3 = c - d*x1_3^2 - x2_3; dx3_3 = r*(s*(x1_3 - x0) - x3_3); dphi_3 = k1*x1_3 - k2*phi_3;
        x1_4 = x1 + dt*dx1_3; x2_4 = x2 + dt*dx2_3; x3_4 = x3 + dt*dx3_3; phi_4 = phi + dt*dphi_3; rho4 = alpha*phi_4^2 + beta*phi_4 + gamma; dx1_4 = x2_4 - a*x1_4^3 + b*x1_4^2 - x3_4 + I_ext - k*rho4*x1_4; dx2_4 = c - d*x1_4^2 - x2_4; dx3_4 = r*(s*(x1_4 - x0) - x3_4); dphi_4 = k1*x1_4 - k2*phi_4;
        
        x1 = x1 + (dt/6)*(dx1_1 + 2*dx1_2 + 2*dx1_3 + dx1_4);
        x2 = x2 + (dt/6)*(dx2_1 + 2*dx2_2 + 2*dx2_3 + dx2_4);
        x3 = x3 + (dt/6)*(dx3_1 + 2*dx3_2 + 2*dx3_3 + dx3_4);
        phi = phi + (dt/6)*(dphi_1 + 2*dphi_2 + 2*dphi_3 + dphi_4);
        
        x1_obs(i) = x1;
    end
    
    % 第三阶段：寻找峰值并计算 ISI
    idx_peaks = find(x1_obs(2:end-1) > x1_obs(1:end-2) & x1_obs(2:end-1) > x1_obs(3:end)) + 1;
    valid_peaks_idx = idx_peaks(x1_obs(idx_peaks) > 0);
    
    if length(valid_peaks_idx) > 1
        % 计算峰间距 ISI (Inter-Spike Interval)
        peaks_time = valid_peaks_idx * dt;
        isi = diff(peaks_time);
        
        % 取唯一值（允许轻微舍入误差），减小绘图冗余数据量
        isi = unique(round(isi*10)/10); 
        
        bif_I = [bif_I; repmat(I_ext, length(isi), 1)];
        bif_ISI = [bif_ISI; isi'];
    end
    
    % 更新进度条
    if mod(idx, 5) == 0
        waitbar(idx/length(I_range), h, sprintf('进度: %d / %d', idx, length(I_range)));
    end
end % <--- 这个 end 就是第 28 行大 for 循环对应的闭合点

close(h);
fprintf('计算完成！正在生成图像...\n');

% 5. 绘制ISI分岔图 (这部分可能您之前复制时没带上)
figure('Position', [200, 200, 700, 500]);
plot(bif_I, bif_ISI, 'k.', 'MarkerSize', 3);
xlabel('外部电流 I_{ext} (mA)', 'FontWeight', 'bold');
ylabel('ISI', 'FontWeight', 'bold');

set(gca, 'FontSize', 12);
xlim([1 4.5]);
ylim([0 200]);