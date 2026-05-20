% =========================================================================
% 用于毕业设计第2.3节: 四阶Runge-Kutta (RK4) 算法基础验证程序
% 求解经典 Lorenz 混沌系统并绘制三维相图
% =========================================================================
clear; clc; close all;

% 1. 系统参数定义 (Lorenz经典混沌参数)
sigma = 10; 
rho = 28; 
beta = 8/3;

% 2. 积分时间与步长设置
dt = 0.01;          % RK4积分步长
t_end = 50;         % 总仿真时间
N = round(t_end/dt);% 总迭代步数

% 3. 初始化状态变量数组
X = zeros(3, N+1);
X(:,1) = [0.1; 0.1; 0.1]; % 赋予初始条件 (x0, y0, z0)

% 4. 核心: 四阶Runge-Kutta (RK4) 迭代求解过程
for i = 1:N
    % 取当前状态
    x_curr = X(1, i);
    y_curr = X(2, i);
    z_curr = X(3, i);
    
    % 计算 k1 (起点斜率)
    k1_x = sigma * (y_curr - x_curr);
    k1_y = x_curr * (rho - z_curr) - y_curr;
    k1_z = x_curr * y_curr - beta * z_curr;
    
    % 计算 k2 (中点斜率1)
    k2_x = sigma * ((y_curr + 0.5*dt*k1_y) - (x_curr + 0.5*dt*k1_x));
    k2_y = (x_curr + 0.5*dt*k1_x) * (rho - (z_curr + 0.5*dt*k1_z)) - (y_curr + 0.5*dt*k1_y);
    k2_z = (x_curr + 0.5*dt*k1_x) * (y_curr + 0.5*dt*k1_y) - beta * (z_curr + 0.5*dt*k1_z);
    
    % 计算 k3 (中点斜率2)
    k3_x = sigma * ((y_curr + 0.5*dt*k2_y) - (x_curr + 0.5*dt*k2_x));
    k3_y = (x_curr + 0.5*dt*k2_x) * (rho - (z_curr + 0.5*dt*k2_z)) - (y_curr + 0.5*dt*k2_y);
    k3_z = (x_curr + 0.5*dt*k2_x) * (y_curr + 0.5*dt*k2_y) - beta * (z_curr + 0.5*dt*k2_z);
    
    % 计算 k4 (终点斜率)
    k4_x = sigma * ((y_curr + dt*k3_y) - (x_curr + dt*k3_x));
    k4_y = (x_curr + dt*k3_x) * (rho - (z_curr + dt*k3_z)) - (y_curr + dt*k3_y);
    k4_z = (x_curr + dt*k3_x) * (y_curr + dt*k3_y) - beta * (z_curr + dt*k3_z);
    
    % 根据加权平均更新下一时刻状态
    X(1, i+1) = x_curr + (dt/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x);
    X(2, i+1) = y_curr + (dt/6) * (k1_y + 2*k2_y + 2*k3_y + k4_y);
    X(3, i+1) = z_curr + (dt/6) * (k1_z + 2*k2_z + 2*k3_z + k4_z);
end

% 5. 绘制求解出的混沌吸引子相图
figure('Position', [300, 200, 700, 500]);
plot3(X(1,:), X(2,:), X(3,:), 'k', 'LineWidth', 0.8);
grid on;
xlabel('x', 'FontWeight', 'bold');
ylabel('y', 'FontWeight', 'bold');
zlabel('z', 'FontWeight', 'bold');

view(45, 20); % 调整最佳三维观测视角