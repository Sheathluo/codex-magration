%% 星链双站无源雷达（Starlink-PBR）链路预算与回波信噪比仿真
% 基于 Starlink Ku 频段下行参数与地面双站接收配置
clear; clc; close all;

%% 1. 核心物理与系统仿真参数配置（与仿真参数表严格一致）
c = 3e8;                        % 光速 (m/s)
k_B = 1.380649e-23;             % 玻尔兹曼常数 (J/K)
T0 = 290;                       % 标准系统噪声温度 (K)

% 星载发射端参数 (Starlink V1.5 / V2 Mini Ku-band)
fc = 11.325e9;                  % 载波中心频率: 11.325 GHz
lambda = c / fc;                % 波长 (m)
B = 240e6;                      % 单信道有效带宽: 240 MHz
EIRP_dBW = 38;                  % 单波束等效全向辐射功率: 38 dBW
EIRP = 10^(EIRP_dBW / 10);      % 线性值 (W)

% 双站空间几何与目标参数
R_T = 500e3;                    % 卫星至目标斜距: 700 km (对应约50°仰角)
R_R = linspace(1e3, 60e3, 500); % 目标至地面接收机斜距: 1 km ~ 60 km
sigma_b_list = [0.05, 1, 5];    % 目标双站 RCS (m^2): 0.05(无人机), 1(小型机), 5(民航机)

% 地面接收端参数
G_rx_dBi = 33;                  % 监控通道天线增益: 33 dBi
G_rx = 10^(G_rx_dBi / 10);      % 线性值
Fn_dB = 2.0;                    % 接收机噪声系数: 2.0 dB
Fn = 10^(Fn_dB / 10);           % 线性值
L_sys_dB = 3.0;                 % 综合损耗 (大气雨衰 1.5 dB + 失配/波束扫描损耗 1.5 dB)
L_sys = 10^(L_sys_dB / 10);     % 线性值

% 相干积累时间与检测门限
T_int_list = [5e-3, 20e-3, 50e-3, 100e-3]; % 相干积累时间: 5, 20, 50, 100 ms
SNR_thresh_dB = 13.0;           % 典型 CFAR 恒虚警检测门限 (Pd=0.9, Pfa=1e-6)

%% 2. 图1：固定目标 RCS (sigma_b = 5 m^2)，不同相干积累时间下的 SNR-斜距曲线
sigma_ref = 5.0; % 基准目标 RCS = 5 m^2
figure('Color', 'w', 'Position', [100, 100, 720, 480]);
hold on; grid on; box on;

colors = lines(length(T_int_list));
for idx = 1:length(T_int_list)
    T_int = T_int_list(idx);
    
    % 双站雷达相干积累后回波信噪比公式:
    % SNR = (EIRP * G_rx * lambda^2 * sigma_b * T_int) / ((4*pi)^3 * R_T^2 * R_R.^2 * k_B * T0 * Fn * L_sys)
    % 注：此处匹配滤波带宽与噪声带宽均为 B，处理增益 G_p = B * T_int，两者在分子分母中 B 约去，保留 T_int
    numerator = EIRP * G_rx * (lambda^2) * sigma_ref * T_int;
    denominator = ((4 * pi)^3) * (R_T^2) * (R_R.^2) * k_B * T0 * Fn * L_sys;
    SNR_linear = numerator ./ denominator;
    SNR_dB = 10 * log10(SNR_linear);
    
    plot(R_R / 1e3, SNR_dB, 'LineWidth', 2.0, 'Color', colors(idx, :), ...
        'DisplayName', sprintf('T_{int} = %d ms (G_p = %.1f dB)', T_int*1e3, 10*log10(B*T_int)));
end

% 绘制检测门限参考线
yline(SNR_thresh_dB, '--r', 'LineWidth', 1.5, 'DisplayName', 'CFAR 检测门限 (13 dB)');

xlabel('目标至地面接收机斜距 R_R (km)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('相干积累后输出信噪比 SNR (dB)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('星链双站无源雷达回波 SNR 性能曲线 (R_T = %d km, \\sigma_b = %.1f m^2, EIRP = %d dBW)', ...
    R_T/1e3, sigma_ref, EIRP_dBW), 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 9.5);
xlim([1, 60]);
ylim([-5, 45]);

%% 3. 图2：不同目标 RCS 下，最大可探测斜距 R_R,max 随相干积累时间 T_int 的变化
T_int_sweep = linspace(1e-3, 150e-3, 300); % 1 ms ~ 150 ms
SNR_thresh_lin = 10^(SNR_thresh_dB / 10);

figure('Color', 'w', 'Position', [850, 100, 720, 480]);
hold on; grid on; box on;

rcs_labels = {'\sigma_b = 0.05 m^2 (小型无人机)', '\sigma_b = 1.0 m^2 (通用飞行器)', '\sigma_b = 5.0 m^2 (民航客机)'};
line_styles = {'-.', '--', '-'};

for s_idx = 1:length(sigma_b_list)
    sigma_val = sigma_b_list(s_idx);
    
    % 反解满足 SNR >= SNR_thresh 的最大接收机斜距 R_R_max
    R_R_max = sqrt((EIRP * G_rx * (lambda^2) * sigma_val .* T_int_sweep) ./ ...
        (((4 * pi)^3) * (R_T^2) * k_B * T0 * Fn * L_sys * SNR_thresh_lin));
    
    plot(T_int_sweep * 1e3, R_R_max / 1e3, 'LineWidth', 2.2, ...
        'LineStyle', line_styles{s_idx}, 'DisplayName', rcs_labels{s_idx});
end

xlabel('相干积累时间 T_{int} (ms)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('最大可探测斜距 R_{R,max} (km)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('不同目标 RCS 下的最大探测距离边界 (SNR_{min} = %.1f dB, R_T = %d km)', ...
    SNR_thresh_dB, R_T/1e3), 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 10);
xlim([1, 150]);

%% 4. 控制台输出关键工作点验证数值
R_check = 10e3;   % 10 km
T_check = 20e-3;  % 20 ms
SNR_check = 10 * log10((EIRP * G_rx * (lambda^2) * sigma_ref * T_check) / ...
    (((4 * pi)^3) * (R_T^2) * (R_check^2) * k_B * T0 * Fn * L_sys));
fprintf('=== 典型工作点链路预算校核 ===\n');
fprintf('卫星斜距 R_T = %.0f km, 目标斜距 R_R = %.0f km, RCS = %.1f m^2, 积累时间 = %.0f ms\n', ...
    R_T/1e3, R_check/1e3, sigma_ref, T_check*1e3);
fprintf('输出相干信噪比 SNR = %.2f dB\n', SNR_check);