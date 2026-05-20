% =========================================================================
% 第三章：修正后的改进四变量Hindmarsh-Rose神经元模型 RK4数值模拟
% 说明：摒弃视觉修饰，完全依靠改变注入能量参数，真实呈现三种不同动力学状态
% =========================================================================
clear; clc; close all;

% 全局参数设置 (标准H-R模型参数)
a = 1; b = 3; c = 1; d = 5; s = 4; r = 0.006; x0 = -1.6;
dt = 0.01;      % 积分步长
t_total = 1000; % 总仿真时间
N = round(t_total/dt); % 总迭代步数
t = (0:N-1)*dt;

%% 1. 非周期电流作用下的动力学 (呈现三种不同状态)
% 状态1: 1.5 mA (簇放电 Bursting)
% 状态2: 3.5 mA (连续高频放电 Tonic Spiking)
% 状态3: 6.0 mA (振荡死亡/静息态 Oscillation Death)
k = 0.4; k1 = 1; k2 = 0.5; alpha = 0.1; beta = 0.02; gamma = 0.2;
I_ext_list = [1.5, 3.5, 6.0]; % mA

figure('Name', '非周期电流下的膜电位时间序列', 'Position', [100 100 800 800]);
for idx = 1:length(I_ext_list)
    I_ext = I_ext_list(idx);
    
    X = zeros(4, N);
    X(:,1) = [-1.5; 0; 0; 0]; % 初始条件
    
    for i = 1:N-1
        k_rk1 = HR_derivatives(X(:,i), I_ext, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0);
        k_rk2 = HR_derivatives(X(:,i) + 0.5*dt*k_rk1, I_ext, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0);
        k_rk3 = HR_derivatives(X(:,i) + 0.5*dt*k_rk2, I_ext, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0);
        k_rk4 = HR_derivatives(X(:,i) + dt*k_rk3, I_ext, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0);
        X(:,i+1) = X(:,i) + (dt/6) * (k_rk1 + 2*k_rk2 + 2*k_rk3 + k_rk4);
    end
    
    subplot(3,1,idx);
    plot(t, X(1,:), 'k', 'LineWidth', 1);
    title(['外部电流 I_{ext} = ', num2str(I_ext), ' mA']);
    xlabel('时间'); ylabel('膜电位 (mV)');
    xlim([0 1000]);
end

%% 2. 周期电流作用下的动力学 (呈现三种不同状态)
% 状态1: A = 0.5 (亚阈值微弱响应 Subthreshold)
% 状态2: A = 2.0 (周期性簇放电 Periodic Bursting)
% 状态3: A = 5.0 (连续密集放电 Tonic Spiking)
k_p = 1; k1_p = 0.9; k2_p = 0.5; alpha_p = 0.4; beta_p = 0.02; gamma_p = 0.2;
omega = 0.02;
A_list = [0.5, 2.0, 5.0];

figure('Name', '周期电流下的膜电位时间序列', 'Position', [150 150 800 800]);
for idx = 1:length(A_list)
    A = A_list(idx);
    
    X = zeros(4, N);
    X(:,1) = [-1.5; 0; 0; 0];
    
    for i = 1:N-1
        I_t1 = A * cos(omega * t(i));
        I_t2 = A * cos(omega * (t(i) + 0.5*dt));
        I_t3 = A * cos(omega * (t(i) + dt));
        
        k_rk1 = HR_derivatives(X(:,i), I_t1, k_p, k1_p, k2_p, alpha_p, beta_p, gamma_p, a, b, c, d, r, s, x0);
        k_rk2 = HR_derivatives(X(:,i) + 0.5*dt*k_rk1, I_t2, k_p, k1_p, k2_p, alpha_p, beta_p, gamma_p, a, b, c, d, r, s, x0);
        k_rk3 = HR_derivatives(X(:,i) + 0.5*dt*k_rk2, I_t2, k_p, k1_p, k2_p, alpha_p, beta_p, gamma_p, a, b, c, d, r, s, x0);
        k_rk4 = HR_derivatives(X(:,i) + dt*k_rk3, I_t3, k_p, k1_p, k2_p, alpha_p, beta_p, gamma_p, a, b, c, d, r, s, x0);
        X(:,i+1) = X(:,i) + (dt/6) * (k_rk1 + 2*k_rk2 + 2*k_rk3 + k_rk4);
    end
    
    subplot(3,1,idx);
    plot(t, X(1,:), 'k', 'LineWidth', 1);
    title(['A = ', num2str(A), ', \omega = 0.02']);
    xlabel('时间'); ylabel('膜电位 (mV)');
    xlim([0 1000]);
end

%% 3. 噪声作用下的动力学
D0 = 0.1;
I_ext_noise = 4.2;

figure('Name', '噪声+非周期电流', 'Position', [200 200 600 300]);
X = zeros(4, N);
X(:,1) = [-1.5; 0; 0; 0];

for i = 1:N-1
    k_rk1 = HR_derivatives(X(:,i), I_ext_noise, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0);
    k_rk2 = HR_derivatives(X(:,i) + 0.5*dt*k_rk1, I_ext_noise, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0);
    k_rk3 = HR_derivatives(X(:,i) + 0.5*dt*k_rk2, I_ext_noise, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0);
    k_rk4 = HR_derivatives(X(:,i) + dt*k_rk3, I_ext_noise, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0);
    
    dX = (dt/6) * (k_rk1 + 2*k_rk2 + 2*k_rk3 + k_rk4);
    X(:,i+1) = X(:,i) + dX;
    X(4,i+1) = X(4,i+1) + sqrt(2 * D0 * dt) * randn;
end

plot(t, X(1,:), 'k', 'LineWidth', 1);
title(['外部电流 I_{ext} = 4.2 mA, D_0 = ', num2str(D0)]);
xlabel('时间'); ylabel('膜电位 (mV)');
xlim([0 1000]);

%% 内部调用的导数计算函数
function dX = HR_derivatives(X, I_ext, k, k1, k2, alpha, beta, gamma, a, b, c, d, r, s, x0)
    x1 = X(1); x2 = X(2); x3 = X(3); phi = X(4);
    rho_phi = alpha * phi^2 + beta * phi + gamma;
    
    dx1 = x2 - a*x1^3 + b*x1^2 - x3 + I_ext - k * rho_phi * x1;
    dx2 = c - d*x1^2 - x2;
    dx3 = r*(s*(x1 - x0) - x3);
    dphi = k1*x1 - k2*phi;
    
    dX = [dx1; dx2; dx3; dphi];
end