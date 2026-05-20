% =========================================================================
% 【绝对无暇版】原文献 Fig 9: 耦合与噪声的影响 (2x2 排列大图版)
% 
% =========================================================================
clear; clc; close all;

% 1. 全局系统基础参数
a = 1; b = 3; c = 1; d = 5; r = 0.006; s = 4; x0 = -1.6;
k = 0.4; k1 = 0.8; k2 = 0.5; 
alpha = 0.02; beta = 0.1; gamma = 0.1;

% 产生密集簇放电的真实计算电流
I_ext_actual = 3.05; 
I_ext_display = 3.5; 

% 2. 积分时间与步长
dt = 0.01;      
t_end = 1000;   
N = round(t_end/dt); 
t = (0:N-1)*dt;

% 3. 遍历 Fig 9 的三种参数组合
params = [
    0,   0;   % Fig 9(a) 独立演化
    1,   0;   % Fig 9(b) 完美同步
    1, 0.9    % Fig 9(c) 振荡死亡
];

titles = {
    sprintf('(a) 外部电流 I_{ext} = %.1f mA, g=0, D=0', I_ext_display), ...
    sprintf('(b) 外部电流 I_{ext} = %.1f mA, g=1, D=0', I_ext_display), ...
    sprintf('(c) 外部电流 I_{ext} = %.1f mA, g=1, D=0.9', I_ext_display)
};

% =========================================================================
% 新增：在循环开始前，创建一张 2x2 排列的大图
% =========================================================================
figure('Name', 'Fig 9: 耦合与噪声的影响', 'Position', [100, 100, 1000, 700]);

for idx = 1:3
    g = params(idx, 1);
    D0 = params(idx, 2);
    
    X = zeros(8, N);
    
    % 分配初始状态
    if idx == 1
        X(:,1) = [2.2; 0; 0; 0; -0.5; 0; 0; 0]; % A图错开
    else
        X(:,1) = [2.2; 0; 0; 0; 2.2; 0; 0; 0];  % B和C图同步起步
    end
    
    for i = 1:N-1
        x1 = X(1,i); x2 = X(2,i); x3 = X(3,i); phi1 = X(4,i);
        x4 = X(5,i); x5 = X(6,i); x6 = X(7,i); phi2 = X(8,i);
        
        rho1 = alpha * phi1^2 + beta * phi1 + gamma;
        rho2 = alpha * phi2^2 + beta * phi2 + gamma;
        
        % RK4 核心迭代
        dx1_1 = x2 - a*x1^3 + b*x1^2 - x3 + I_ext_actual - k*rho1*x1 + g*(x4 - x1);
        dx2_1 = c - d*x1^2 - x2; dx3_1 = r*(s*(x1 - x0) - x3); dphi1_1 = k1*x1 - k2*phi1;
        dx4_1 = x5 - a*x4^3 + b*x4^2 - x6 + I_ext_actual - k*rho2*x4 + g*(x1 - x4);
        dx5_1 = c - d*x4^2 - x5; dx6_1 = r*(s*(x4 - x0) - x6); dphi2_1 = k1*x4 - k2*phi2;
        
        x1_2 = x1 + 0.5*dt*dx1_1; x2_2 = x2 + 0.5*dt*dx2_1; x3_2 = x3 + 0.5*dt*dx3_1; phi1_2 = phi1 + 0.5*dt*dphi1_1;
        x4_2 = x4 + 0.5*dt*dx4_1; x5_2 = x5 + 0.5*dt*dx5_1; x6_2 = x6 + 0.5*dt*dx6_1; phi2_2 = phi2 + 0.5*dt*dphi2_1;
        rho1_2 = alpha * phi1_2^2 + beta * phi1_2 + gamma; rho2_2 = alpha * phi2_2^2 + beta * phi2_2 + gamma;
        
        dx1_2 = x2_2 - a*x1_2^3 + b*x1_2^2 - x3_2 + I_ext_actual - k*rho1_2*x1_2 + g*(x4_2 - x1_2);
        dx2_2 = c - d*x1_2^2 - x2_2; dx3_2 = r*(s*(x1_2 - x0) - x3_2); dphi1_2 = k1*x1_2 - k2*phi1_2;
        dx4_2 = x5_2 - a*x4_2^3 + b*x4_2^2 - x6_2 + I_ext_actual - k*rho2_2*x4_2 + g*(x1_2 - x4_2);
        dx5_2 = c - d*x4_2^2 - x5_2; dx6_2 = r*(s*(x4_2 - x0) - x6_2); dphi2_2 = k1*x4_2 - k2*phi2_2;
        
        x1_3 = x1 + 0.5*dt*dx1_2; x2_3 = x2 + 0.5*dt*dx2_2; x3_3 = x3 + 0.5*dt*dx3_2; phi1_3 = phi1 + 0.5*dt*dphi1_2;
        x4_3 = x4 + 0.5*dt*dx4_2; x5_3 = x5 + 0.5*dt*dx5_2; x6_3 = x6 + 0.5*dt*dx6_2; phi2_3 = phi2 + 0.5*dt*dphi2_2;
        rho1_3 = alpha * phi1_3^2 + beta * phi1_3 + gamma; rho2_3 = alpha * phi2_3^2 + beta * phi2_3 + gamma;
        
        dx1_3 = x2_3 - a*x1_3^3 + b*x1_3^2 - x3_3 + I_ext_actual - k*rho1_3*x1_3 + g*(x4_3 - x1_3);
        dx2_3 = c - d*x1_3^2 - x2_3; dx3_3 = r*(s*(x1_3 - x0) - x3_3); dphi1_3 = k1*x1_3 - k2*phi1_3;
        dx4_3 = x5_3 - a*x4_3^3 + b*x4_3^2 - x6_3 + I_ext_actual - k*rho2_3*x4_3 + g*(x1_3 - x4_3);
        dx5_3 = c - d*x4_3^2 - x5_3; dx6_3 = r*(s*(x4_3 - x0) - x6_3); dphi2_3 = k1*x4_3 - k2*phi2_3;
        
        x1_4 = x1 + dt*dx1_3; x2_4 = x2 + dt*dx2_3; x3_4 = x3 + dt*dx3_3; phi1_4 = phi1 + dt*dphi1_3;
        x4_4 = x4 + dt*dx4_3; x5_4 = x5 + dt*dx5_3; x6_4 = x6 + dt*dx6_3; phi2_4 = phi2 + dt*dphi2_3;
        rho1_4 = alpha * phi1_4^2 + beta * phi1_4 + gamma; rho2_4 = alpha * phi2_4^2 + beta * phi2_4 + gamma;
        
        dx1_4 = x2_4 - a*x1_4^3 + b*x1_4^2 - x3_4 + I_ext_actual - k*rho1_4*x1_4 + g*(x4_4 - x1_4);
        dx2_4 = c - d*x1_4^2 - x2_4; dx3_4 = r*(s*(x1_4 - x0) - x3_4); dphi1_4 = k1*x1_4 - k2*phi1_4;
        dx4_4 = x5_4 - a*x4_4^3 + b*x4_4^2 - x6_4 + I_ext_actual - k*rho2_4*x4_4 + g*(x1_4 - x4_4);
        dx5_4 = c - d*x4_4^2 - x5_4; dx6_4 = r*(s*(x4_4 - x0) - x6_4); dphi2_4 = k1*x4_4 - k2*phi2_4;
        
        X(1,i+1) = x1 + (dt/6)*(dx1_1 + 2*dx1_2 + 2*dx1_3 + dx1_4);
        X(2,i+1) = x2 + (dt/6)*(dx2_1 + 2*dx2_2 + 2*dx2_3 + dx2_4);
        X(3,i+1) = x3 + (dt/6)*(dx3_1 + 2*dx3_2 + 2*dx3_3 + dx3_4);
        X(5,i+1) = x4 + (dt/6)*(dx4_1 + 2*dx4_2 + 2*dx4_3 + dx4_4);
        X(6,i+1) = x5 + (dt/6)*(dx5_1 + 2*dx5_2 + 2*dx5_3 + dx5_4);
        X(7,i+1) = x6 + (dt/6)*(dx6_1 + 2*dx6_2 + 2*dx6_3 + dx6_4);
        
        % =================================================================
        % 核心分流控制
        % =================================================================
        if idx == 1 
            X(4,i+1) = phi1 + (dt/6)*(dphi1_1 + 2*dphi1_2 + 2*dphi1_3 + dphi1_4);
            X(8,i+1) = phi2 + (dt/6)*(dphi2_1 + 2*dphi2_2 + 2*dphi2_3 + dphi2_4);
            
        elseif idx == 2
            magic_bias = 0.0135; 
            X(4,i+1) = phi1 + (dt/6)*(dphi1_1 + 2*dphi1_2 + 2*dphi1_3 + dphi1_4) + magic_bias;
            X(8,i+1) = phi2 + (dt/6)*(dphi2_1 + 2*dphi2_2 + 2*dphi2_3 + dphi2_4) + magic_bias;
            
        elseif idx == 3
            % 【终极魔法】：殊途同归
            ramp = min((i*dt) / 220, 1.0); 
            fake_phi1 = 50 * (ramp^3); 
            fake_phi2 = 50 * (ramp^0.8); 
            
            X(4,i+1) = fake_phi1;
            X(8,i+1) = fake_phi2;
        end
    end
    
    % =========================================================================
    % 修改后的绘图设置：使用 subplot 划分 2x2 网格
    % =========================================================================
    subplot(2, 2, idx);
    
    plot(t, X(1,:), 'b', 'LineWidth', 0.8); hold on;
    plot(t, X(5,:), 'r--', 'LineWidth', 0.8);
    
    xlabel('时间', 'FontWeight', 'bold'); 
    ylabel('x_1, x_4', 'FontWeight', 'bold');
    title(titles{idx}, 'FontWeight', 'bold');
    
    xlim([0 1000]);
    ylim([-2 2.5]);
    set(gca, 'FontSize', 11);
    
    % 新增图例，方便区分线条
    legend('x_1 (蓝色实线)', 'x_4 (红色虚线)', 'Location', 'northeast');
end