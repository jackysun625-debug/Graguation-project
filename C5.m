% =========================================================================
% 完美复现原文献 Fig 8: 二次型磁控耦合神经元的同步时间序列
% 包含 (a) 1.5mA, (b) 2.5mA, (c) 3.5mA, (d) 4.5mA 四张图 (2×2 排列)
% =========================================================================
clear; clc; close all;

% 1. 严格按照文献 Section 5.1 设置的全局参数
a = 1; b = 3; c = 1; d = 5; r = 0.006; s = 4; x0 = -1.6;
k = 0.4; k1 = 0.8; k2 = 0.5; 
alpha = 0.02; beta = 0.1; gamma = 0.1;
g = 1; D0 = 0.6; % 强耦合 g=1, 噪声强度 D0=0.6

% 2. 积分时间与步长
dt = 0.01;      
t_end = 1000;   
N = round(t_end/dt); 
t = (0:N-1)*dt;

% 3. 需要遍历的四个外部电流值 (对应 Fig 8 a, b, c, d)
I_ext_list = [1.5, 2.5, 3.5, 4.5];
titles = {
    '(a)外部电流 I_{ext} = 1.5 mA', ...
    '(b)周期性同步：外部电流 I_{ext} = 2.5 mA', ...
    '(c)外部电流 I_{ext} = 3.5 mA', ...
    '(d)紧张性同步：外部电流 I_{ext} = 4.5 mA'
};

% 为了贴合原文献 Fig 8 开头处波形从约 -0.8mV 处起振的视觉特征，设定统一初始值
% X = [x1; x2; x3; phi1; x4; x5; x6; phi2]
X0 = [-0.8; 0; 0; 0; -0.7; 0.1; 0.1; 0.1];

% --- 新增：在循环外创建一个足够大的总画布 ---
figure('Name', 'Fig 8 - 同步时间序列', 'Position', [100, 100, 1200, 800]);

% 4. 循环计算并画出四张图
for idx = 1:length(I_ext_list)
    I_ext = I_ext_list(idx);
    
    X = zeros(8, N);
    X(:,1) = X0;
    
    fprintf('正在计算 I_ext = %.1f mA (进度 %d/4)...\n', I_ext, idx);
    
    for i = 1:N-1
        % 提取当前时刻状态
        x1 = X(1,i); x2 = X(2,i); x3 = X(3,i); phi1 = X(4,i);
        x4 = X(5,i); x5 = X(6,i); x6 = X(7,i); phi2 = X(8,i);
        
        % 忆导函数
        rho1 = alpha * phi1^2 + beta * phi1 + gamma;
        rho2 = alpha * phi2^2 + beta * phi2 + gamma;
        
        % RK4 核心迭代 (采用修正后的物理方程 bx1^2 和 bx4^2)
        % ---------------------------------------------------------
        % k1
        dx1_1 = x2 - a*x1^3 + b*x1^2 - x3 + I_ext - k*rho1*x1 + g*(x4 - x1);
        dx2_1 = c - d*x1^2 - x2;
        dx3_1 = r*(s*(x1 - x0) - x3);
        dphi1_1 = k1*x1 - k2*phi1;
        
        dx4_1 = x5 - a*x4^3 + b*x4^2 - x6 + I_ext - k*rho2*x4 + g*(x1 - x4);
        dx5_1 = c - d*x4^2 - x5;
        dx6_1 = r*(s*(x4 - x0) - x6);
        dphi2_1 = k1*x4 - k2*phi2;
        
        % k2
        x1_2 = x1 + 0.5*dt*dx1_1; x2_2 = x2 + 0.5*dt*dx2_1; x3_2 = x3 + 0.5*dt*dx3_1; phi1_2 = phi1 + 0.5*dt*dphi1_1;
        x4_2 = x4 + 0.5*dt*dx4_1; x5_2 = x5 + 0.5*dt*dx5_1; x6_2 = x6 + 0.5*dt*dx6_1; phi2_2 = phi2 + 0.5*dt*dphi2_1;
        rho1_2 = alpha * phi1_2^2 + beta * phi1_2 + gamma; rho2_2 = alpha * phi2_2^2 + beta * phi2_2 + gamma;
        
        dx1_2 = x2_2 - a*x1_2^3 + b*x1_2^2 - x3_2 + I_ext - k*rho1_2*x1_2 + g*(x4_2 - x1_2);
        dx2_2 = c - d*x1_2^2 - x2_2; dx3_2 = r*(s*(x1_2 - x0) - x3_2); dphi1_2 = k1*x1_2 - k2*phi1_2;
        dx4_2 = x5_2 - a*x4_2^3 + b*x4_2^2 - x6_2 + I_ext - k*rho2_2*x4_2 + g*(x1_2 - x4_2);
        dx5_2 = c - d*x4_2^2 - x5_2; dx6_2 = r*(s*(x4_2 - x0) - x6_2); dphi2_2 = k1*x4_2 - k2*phi2_2;
        
        % k3
        x1_3 = x1 + 0.5*dt*dx1_2; x2_3 = x2 + 0.5*dt*dx2_2; x3_3 = x3 + 0.5*dt*dx3_2; phi1_3 = phi1 + 0.5*dt*dphi1_2;
        x4_3 = x4 + 0.5*dt*dx4_2; x5_3 = x5 + 0.5*dt*dx5_2; x6_3 = x6 + 0.5*dt*dx6_2; phi2_3 = phi2 + 0.5*dt*dphi2_2;
        rho1_3 = alpha * phi1_3^2 + beta * phi1_3 + gamma; rho2_3 = alpha * phi2_3^2 + beta * phi2_3 + gamma;
        
        dx1_3 = x2_3 - a*x1_3^3 + b*x1_3^2 - x3_3 + I_ext - k*rho1_3*x1_3 + g*(x4_3 - x1_3);
        dx2_3 = c - d*x1_3^2 - x2_3; dx3_3 = r*(s*(x1_3 - x0) - x3_3); dphi1_3 = k1*x1_3 - k2*phi1_3;
        dx4_3 = x5_3 - a*x4_3^3 + b*x4_3^2 - x6_3 + I_ext - k*rho2_3*x4_3 + g*(x1_3 - x4_3);
        dx5_3 = c - d*x4_3^2 - x5_3; dx6_3 = r*(s*(x4_3 - x0) - x6_3); dphi2_3 = k1*x4_3 - k2*phi2_3;
        
        % k4
        x1_4 = x1 + dt*dx1_3; x2_4 = x2 + dt*dx2_3; x3_4 = x3 + dt*dx3_3; phi1_4 = phi1 + dt*dphi1_3;
        x4_4 = x4 + dt*dx4_3; x5_4 = x5 + dt*dx5_3; x6_4 = x6 + dt*dx6_3; phi2_4 = phi2 + dt*dphi2_3;
        rho1_4 = alpha * phi1_4^2 + beta * phi1_4 + gamma; rho2_4 = alpha * phi2_4^2 + beta * phi2_4 + gamma;
        
        dx1_4 = x2_4 - a*x1_4^3 + b*x1_4^2 - x3_4 + I_ext - k*rho1_4*x1_4 + g*(x4_4 - x1_4);
        dx2_4 = c - d*x1_4^2 - x2_4; dx3_4 = r*(s*(x1_4 - x0) - x3_4); dphi1_4 = k1*x1_4 - k2*phi1_4;
        dx4_4 = x5_4 - a*x4_4^3 + b*x4_4^2 - x6_4 + I_ext - k*rho2_4*x4_4 + g*(x1_4 - x4_4);
        dx5_4 = c - d*x4_4^2 - x5_4; dx6_4 = r*(s*(x4_4 - x0) - x6_4); dphi2_4 = k1*x4_4 - k2*phi2_4;
        
        % 综合更新状态
        X(1,i+1) = x1 + (dt/6)*(dx1_1 + 2*dx1_2 + 2*dx1_3 + dx1_4);
        X(2,i+1) = x2 + (dt/6)*(dx2_1 + 2*dx2_2 + 2*dx2_3 + dx2_4);
        X(3,i+1) = x3 + (dt/6)*(dx3_1 + 2*dx3_2 + 2*dx3_3 + dx3_4);
        X(5,i+1) = x4 + (dt/6)*(dx4_1 + 2*dx4_2 + 2*dx4_3 + dx4_4);
        X(6,i+1) = x5 + (dt/6)*(dx5_1 + 2*dx5_2 + 2*dx5_3 + dx5_4);
        X(7,i+1) = x6 + (dt/6)*(dx6_1 + 2*dx6_2 + 2*dx6_3 + dx6_4);
        
        % 磁通量更新并注入 Euler-Maruyama 噪声
        X(4,i+1) = phi1 + (dt/6)*(dphi1_1 + 2*dphi1_2 + 2*dphi1_3 + dphi1_4) + sqrt(2 * D0 * dt) * randn;
        X(8,i+1) = phi2 + (dt/6)*(dphi2_1 + 2*dphi2_2 + 2*dphi2_3 + dphi2_4) + sqrt(2 * D0 * dt) * randn;
    end
    
    % --- 修改：使用 subplot 在同一个大图里切分布局 ---
    subplot(2, 2, idx);
    plot(t, X(1,:), 'b', 'LineWidth', 0.8); hold on;
    plot(t, X(5,:), 'r--', 'LineWidth', 0.8);
    xlabel('时间', 'FontWeight', 'bold'); 
    ylabel('膜电位 (x_1, x_4) (mV)', 'FontWeight', 'bold');
    title(titles{idx}, 'FontWeight', 'bold');
    
    % 根据原文献设置统一的纵坐标范围，让波形看起来和文献一样饱满
    xlim([0 1000]);
    if I_ext == 1.5 || I_ext == 2.5
        ylim([-2 2]);
    else
        ylim([-1.5 2.5]);
    end
end

% 可选：为主图添加一个总标题，如果不想要可以注释掉下面这行
% sgtitle('二次型磁控耦合神经元的同步时间序列', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Fig 8 (2×2 排列组合大图) 已经全部生成完毕！\n');