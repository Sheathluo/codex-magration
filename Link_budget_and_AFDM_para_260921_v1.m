%%% 时间：260318
% 链路预算
%
%%%%%%%%%%
clear
% clc
% close all

%% 物理常数
c = 3e8;
k = 1.38e-23;                   % 波尔茨曼常数

%% 信号参数
wave_para = 1;
if wave_para == 1
    % 星链下行链路，面向星链专用的相控阵终端
    fc = 11.7e9;
    lambda = c/fc;
    N = 1024;                   % 子载波数
    alpha = 0;               % 滚降因子 3G WCDMA 明确写死：RRC α = 0.22（TS 25.104/25.213 一脉相承）原值0.25
    L = 1;                      % 过采样倍数
    B = 240e6;
    fs = L*B*(1+alpha);
    delta_f = B/N;  % 子载波间隔
    Ncp = 32;
    T_cp = Ncp/fs;
    T_symbol = 1/delta_f + T_cp;       % 符号持续时间
elseif wave_para == 2
    % 星链下行链路，4G标准
    fc = 2e9;
    lambda = c/fc;
    Nfft = 512;             % 4G的fft点数
    N = 300;                % 真正有效子载波数
    delta_f = 15e3;
    B = N*delta_f;
    Ncp = 36;
    alpha = 0;              % 滚降因子
    fs = 7.68e6;            % 512(5 MHz 模式下采样率 7.68MHz) 遵循标准 LTE 采样率 (30.72MHz/4)
    T_cp = Ncp/fs;
    T_symbol = T_cp + 1/delta_f;       % 符号持续时间
else
    fprintf('不支持的波形参数类型');
end

%% 感知目标参数
object = 4;
if object == 1
    % 假设低轨卫星轨道高度350km
    str = '卫星';
    RCS_dBsm = 15; % 低轨卫星RCS
    Range = 235e3; % 低轨卫星感知距离
    velocity_satellit = 7.7e3; % 速度 km/s
    inclination_one = 38;       % 轨道倾角
    inclination_two = 142;      % 轨道倾角
    velocity = sind(inclination_one)*velocity_satellit+...
        velocity_satellit*sind(inclination_two);
    % 环境参数
    Temperature = 393;              % 温度120°，航空级元器件能承受的最高温度
elseif object == 2
    str = '导弹';
    RCS_dBsm = 5; % 弹道导弹RCS
    Range = 200e3; % 弹道导弹飞行高度100km-200km
    velocity = 6.8e3; % 速度 20马赫 km/s
    % 环境参数
    Temperature = 339;              % 温度50～80°，导引头能承受的温度范围
elseif object == 3
    str = '小型无人机';
    RCS_dBsm = -20; % 小型无人机RCS
    Range = 1000; % 小型无人机10m-1km
    velocity = 40; % 相对速度m/s    无人机速度：2km/h ~ 120km/h(33m/s)，平台速度：20km/h
    % 环境参数
    Temperature = 292;              % 室温20度
elseif object == 4
    str = '战斗机';
    RCS_dBsm = 20; % 战斗机RCS
    % Range = 330e3; % 战斗机飞行最高飞行20km
    % 假设卫星高度480km。星链Gen2 V2 Mini ~480km, 星链Gen2手机直连D2C ~350-360km,Gen2
    % 新规划的极低轨道~326-345km
    Range_sate2fighter = 466e3; % 卫星到目标的距离
    Range_fighter2radar = 100e3; % 目标到雷达的距离，一般机载雷达的感知距离
    Range = sqrt(Range_sate2fighter * Range_fighter2radar);
    % Range_seta2rx = sqrt(Range_sate2fighter^2+Range_fighter2radar^2);
    velocity = 0.681e3; % 速度 2马赫
    % 环境参数
    Temperature = 216;              % 高空的温度
    % tau_ = Range / c * 1e3
else
    error("感知目标参数必须是1-4");
end
RCS = 10^(RCS_dBsm/10);         %% 感知目标RCS=10m^2;

P_t_dBm = 42;               % 天线单通道的发射功率，单位dBm
P_t = 10^((P_t_dBm-30)/10); % 发射功率  单位W
%% 天线增益
Gt_dB = 37;                  % 发射端天线带增益(dB)
Gt_linear = 10^(Gt_dB/10);
Gr_dB = 34;                  % 接收端天线带增益(dB)
Gr_linear = 10^(Gr_dB/10);
%% 发射端的EIRP
% EIRP_dBm = P_t_dBm+Gt_dB;                       % 单位dBm EIRP(dBm) = P_T(dBm)+Gt(dB)
EIRP_dBm = 39 + 30;
EIRP_linear = 10^((EIRP_dBm-30)/10);            % 单位W
%% 处理损耗
L_dB = 3;
L = 10^(L_dB/10);
%% 噪声系数
NF_dB = 2;                  % 典型接收机噪声系数1.5~2.5
F = 10^(NF_dB/10);
P_noise = B * k * Temperature * F;                              % 底噪功率
%% 接收端宽带输入SNR（累积处理前）
% P_r = P_t*Gt_linear*Gr_linear*lambda^2*RCS/((4*pi)^3*Range^4*L);      % 输入功率
P_r = EIRP_linear * Gr_linear * lambda^2 * RCS/((4*pi)^3*Range^4*L);    % 输入功率

SNR_wideband_dB = 10*log10(P_r/P_noise);
%% 检测概率为95%，此处的虚警概率按照10^-6设置，理论值为13.9dB
SNR_required_after_accum_dB = 14;               % 累积处理后的目标SNR
G_required_dB = SNR_required_after_accum_dB - SNR_wideband_dB;
G_required = 10^(G_required_dB/10);

% 理想相干累积增益为 B*T；先求最小时间，再向上取整到完整符号
T_required = G_required/B;                        % 理论最小累积时间，单位：s
Nsym = ceil(T_required/T_symbol);                 % 累积符号数量
T_accum = Nsym*T_symbol;                          % 实际累积时间，单位：s
G_actual_dB = 10*log10(B*T_accum);                % 取整后的实际相干处理增益
SNR_after_one_symbol_dB = SNR_wideband_dB + 10*log10(B*T_symbol);
SNR_after_accum_dB = SNR_wideband_dB + G_actual_dB;
%% 雷达系统参数
theta_rad = acos(50/466);           % 双站雷达的双站角度
delta_R = c/(2*B*cos(theta_rad/2));   % 距离分辨率
delta_v = lambda/(2*T_accum*cos(theta_rad/2));           % 速度分辨率
R_bin = Range/delta_R;
v_bin = velocity/delta_v;
N_radar_fig = R_bin*v_bin;

%% 徙动
range_migration = T_accum*velocity;
range_bins = Range/delta_R;
range_migration_bins = range_migration/delta_R;

f_d = 2*velocity*cos(theta_rad/2)/lambda;
norm_doppler = f_d/delta_f;

delay_time = Range_fighter2radar / c;

velocity_bins = velocity/delta_v;
alpha_max = ceil(norm_doppler);
%% 直射径与反射径的相对时延近似计算
theta = 90;
h_r = 60;                       % hr是接收机离地高度 60m 典型无人机的飞行高度
tau = 2*h_r*sind(theta)/c;
norm_tau = tau*fs;
l_max = ceil(norm_tau);

%% 直射径与反射径的相对时延近似计算
% los_path = 466e3;      % 直射径
% reflect_path = 466e3+100e3;
% tau = (reflect_path-los_path)/c;
% norm_tau = tau*fs;
% l_max = ceil(norm_tau);
%% 参数显示
fprintf('=== 系统参数 ===\n');
fprintf('载频：%.1f GHz \n', fc/1e9);
fprintf('带宽：%.3f MHz \n', B/1e6);
fprintf('采样率：%.2f MHz \n', fs/1e6);
fprintf('符号持续时间：%.2f us \n', T_symbol*1e6);

fprintf('=== 目标参数 ===\n');
fprintf('目标类型： %s \n', str);
% fprintf('卫星1的轨道倾角：%d° \n', inclination_one);
% fprintf('卫星2的轨道倾角：%d° \n', inclination_two);
fprintf('目标距离：%.1f km \n', Range/1e3);
fprintf('目标速度：%.2f m/s \n', velocity);
fprintf('RCS：%.3f m^2 \n', RCS);

% fprintf('=== 通信接收端参数 ===\n');
% fprintf('卫星仰角：%.1f° \n', theta')
% fprintf('通信接收端离地表高度：%.1f m \n', h_r);
% fprintf('两径真实路程差：%.2f m \n', tao*c);
% fprintf('两径真实时延差：%.2f ns \n', tao*1e9);
% fprintf('归一化时延：%.4f \n', norm_tao);

fprintf('=== 感知性能 ===\n');
fprintf('距离分辨率：%.3f m \n', delta_R);
fprintf('速度分辨率：%.3f m/s \n', delta_v);

fprintf('=== 累积时间计算 ===\n');
fprintf('累积处理前的宽带输入SNR：%.1f dB \n', SNR_wideband_dB);
fprintf('单符号匹配处理后的SNR：%.1f dB \n', SNR_after_one_symbol_dB);
fprintf('累积处理后的目标SNR：%.1f dB \n', SNR_required_after_accum_dB);
fprintf('相干处理增益需求：%.2f dB \n', G_required_dB);
fprintf('理论最小累积时间：%.3f ms \n', T_required*1e3);
fprintf('累积时间对应的符号数量：%d \n', Nsym);
fprintf('取整后的实际累积时间：%.3f ms \n', T_accum*1e3);
fprintf('取整后的实际累积SNR：%.1f dB \n', SNR_after_accum_dB);
fprintf('累积时间内目标位移：%.1f m \n', range_migration);

fprintf('==== 徙动 ====\n');
fprintf('距离对应的距离单元：%.1f \n', range_bins);
fprintf('累积时间内距离徙动对应的距离单元：%.1f \n', range_migration_bins);
fprintf('多普勒频偏：%.4f MHz \n', f_d/1e6);
fprintf('归一化多普勒：%.4f \n', norm_doppler);
fprintf('速度对应的速度单元：%.1f  \n', velocity_bins);

%% AFDM波形参数
xi_v = 4;                       % 扩散间隔
%% Ncp设计
Ncp_low = ceil(norm_tau);
Ncp = ceil(N/32);
%% N设计
N_low = (2*xi_v+1)*(ceil(norm_tau)+1)*B/(B-2*f_d*(ceil(norm_tau)+1));
%% C1设计
C1_low = f_d/B+(2*xi_v+1)/(2*N);
C1_up = 1/2*(ceil(norm_tau)+1);
c1 = (2*(ceil(norm_doppler)+xi_v)+1)/(2*N);

%% 保护间隔
Q = (l_max+1)*(2*(alpha_max+xi_v)+1)-1;
%% 系统多普勒容限
f_d_max = alpha_max * delta_f;
%% AFDM波形参数设计
fprintf('=== AFDM波形参数 ===\n');
fprintf('N = %d \n', N);
fprintf('子载波间隔 = %.3f kHz \n', delta_f/1e3);
fprintf('c1 = %.5f \n', c1);
fprintf('2Nc1 = %.1f \n', 2*N*c1);
fprintf('xi_v = %d \n', xi_v);
fprintf('N_cp = %d \n', Ncp);     % CP长度设计
fprintf('过采样倍数 = %.1f \n', fs/N/delta_f);
fprintf('滚降因子 = %.3f \n', alpha);
%% AFDM波形参数验证
fprintf('=== AFDM波形参数验证 ===\n');
fprintf('N下界: %.1f \n', N_low);
fprintf('Ncp下界: %.1f \n', Ncp_low);

if c1 >= C1_low && c1 <= C1_up
    fprintf('c1 = %.4f 在区间 [%.4f, %.1f] 内\n', c1, C1_low, C1_up);
else
    fprintf('c1 = %.4f 不在区间 [%.4f, %.1f] 内\n', c1, C1_low, C1_up);
end

if Ncp >= Ncp_low
    fprintf('Ncp = %.d 大于 Ncp的下界 = %.1f \n', Ncp, Ncp_low);
else
    fprintf('Ncp = %.d 小于 Ncp的下界 = %.1f ！！！ \n', Ncp, Ncp_low);
end

if N >= Ncp_low
    fprintf('N = %.d 大于 Ncp的下界 = %.1f \n', N, Ncp_low);
else
    fprintf('N = %.d 小于 Ncp的下界 = %.1f ！！！ \n', N, Ncp_low);
end

if N >= N_low
    fprintf('N = %.d 大于 N的下界 = %.1f \n', N, N_low);
else
    fprintf('N = %.d 小于 N的下界 = %.1f ！！！ \n', N, N_low);
end
