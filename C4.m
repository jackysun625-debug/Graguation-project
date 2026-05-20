% =========================================================================
% 改进H-R模型 哈密顿能量分析绘图 (Average Hamilton Energy vs I_ext)
% 严格匹配文献 Fig. 7b: Variation of average Hamilton energy
% =========================================================================
clear; clc; close all;

% 1. 模型全局参数设置
a = 1; b = 3; c = 1; d = 5; s = 4; r = 0.006; x0 = -1.6;
k = 0.4; k1 = 1; k2 = 0.5; alpha = 0.1; beta = 0.02; gamma = 0.2;

dt = 0.02;
t_transient = 1000; 
t_observe = 1000;
N_trans = round(t_transient/dt);
N_obs = round(t_observe/dt);

I_range = 1.0 : 0.05 : 6.0; % 扫描电流范围
avg_energy = zeros(size(I_range));

fprintf('哈密顿能量计算中...\n');

for idx = 1:length(I_range)
    I_ext = I_range(idx);
    x1 = -1.5; x2 = 0; x3 = 0; phi = 0;
    
    % 跑掉瞬态
    for i = 1:N_trans
        rho = alpha*phi^2 + beta*phi + gamma;
        dx1 = x2 - a*x1^3 + b*x1^2 - x3 + I_ext - k*rho*x1;
        dx2 = c - d*x1^2 - x2;
        dx3 = r*(s*(x1 - x0) - x3);
        dphi = k1*x1 - k2*phi;
        x1 = x1 + dx1*dt; x2 = x2 + dx2*dt; x3 = x3 + dx3*dt; phi = phi + dphi*dt;
    end
    
    % 观察期并计算能量
    E_sum = 0;
    for i = 1:N_obs
        rho = alpha*phi^2 + beta*phi + gamma;
        dx1 = x2 - a*x1^3 + b*x1^2 - x3 + I_ext - k*rho*x1;
        dx2 = c - d*x1^2 - x2;
        dx3 = r*(s*(x1 - x0) - x3);
        dphi = k1*x1 - k2*phi;
        x1 = x1 + dx1*dt; x2 = x2 + dx2*dt; x3 = x3 + dx3*dt; phi = phi + dphi*dt;
        
        % 哈密顿能量公式
        H = (2/3)*d*x1^3 - 2*c*x1 + r*s*(x1 - x0)^2 + (x2 - x3 + I_ext - phi)^2 + k2*x1^2;
        E_sum = E_sum + H;
    end
    avg_energy(idx) = E_sum / N_obs;
end

% 绘图
figure;
plot(I_range, avg_energy, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('外部电流 I_{ext}');
ylabel('哈密顿能量');

grid on;