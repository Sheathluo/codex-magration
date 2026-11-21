%%%%%%%%%%%%%%%%%%%%%%%%%
% 版本信息：main_20261104_v11_Monte_MagriCom_AFDM_N128
% 函数功能：仿真主函数:1.只跑一次，不仿真多次
% 创建日期：2026年11月4日
% 创建人：LuoTao
%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%% ===================== 0. 信号基本参数配置 =====================
para = struct();
% 物理常量
para.c = 3e8;

% 系统基本参数
para.fc = 24e9;
para.lambda = para.c/para.fc;
para.fs_up = 10e6;              % 上采样后采样率
para.L = 2;                     % 过采样倍数
para.L_delay = para.L;          % 理论模糊函数的时延维过采样倍数
para.L_doppler = para.L;        % 理论模糊函数/接收搜索的多普勒维过采样倍数
para.fs = para.fs_up / para.L;  % 采样率
para.Ts = 1/para.fs_up;         % 采样间隔
para.B = para.fs;               % 带宽
para.sinc_filter_len = 300;       % 插值滤波器阶数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para.Nsym = 500;                % 符号数/帧数
add_noise_flag = true;         % 是否加噪 true  false
para.SNR = 5;
% 随机种子
% rng('shuffle');
% rng(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 脉冲成形参数
para.beta = 0.35;               % RRC/RC 滚降因子
% para.beta = 0.25;               % RRC/RC 滚降因子
para.N_preSuffix = 20;          % 每个 AFDM 符号前后额外保护点长度
para.N_tapFiler = 20;           % 成形滤波器半边抽头参数

% AFDM 基本参数
para.Nc = 128;                  % 子载波数
% para.Nc = 1024;                  % 子载波数
para.c1 = 8 / (2 * para.Nc);    % c1
% para.c1 = 11 / (2 * para.Nc);    % c1
para.c2 = 0;                    % c2
para.Mod = 16;                  % 16-QAM
para.mu_4 = 1.32;               % 16QAM四阶矩
para.N = para.Nc;
para.M = para.Nc;
para.delta_f = para.B / para.Nc;    % 子载波间隔
para.T = 1 / para.delta_f;          % 有效符号持续时间
para.Ncp = 20;                      % 循环前缀点数(上采样前)
para.Tcp = para.Ncp / para.B;       % CP持续时间
para.T_AFDM = para.T + para.Tcp;    % AFDM符号持续时间
para.T_symbol = (para.Nc + para.Ncp + 2 * para.N_preSuffix) * para.L / para.fs_up; % 相邻有效符号起点的实际时间间隔

% 加前后保护间隔的脉冲成形后的AFDM符号长度(采样点)
para.samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;
% 发送信号总长度
para.len_signal = para.samples_per_symbol * para.L * para.Nsym + para.N_preSuffix * para.L;

fprintf('\n=============== 信号基本参数 ==============\n');
fprintf('信号带宽 B          = %.1f MHz\n', para.B/1e6);
fprintf('采样率  fs         = %.1f Msps\n', para.fs/1e6);
fprintf('符号数  Nsym       = %d \n', para.Nsym);
fprintf('感知累积时间  T         = %.4f s\n', para.T_symbol*para.Nsym);
fprintf('完整AFDM符号长(采样点)   = %d \n', para.samples_per_symbol);
fprintf('完整信号长度(采样点)      = %d \n', para.len_signal);
fprintf('\n============================================\n');

%% ===================== 1. 感知基本参数配置 ====================
cfg = struct();

% 接收截取与搜索
cfg.index_start = 0;            % 接收信号起始截取点
% cfg.rx_L_delay = 1;           % 接收端时延搜索步进为 1/rx_L_delay

% 感知信道目标信道参数
cfg.Max_num_target_est = 1;       % 接收端 OMP 最多检测目标数，不一定等于真实目标数
cfg.num_target = 1;             % 目标数
cfg.max_nor_doppler_abs = 15;       % 接收端粗搜索最大归一化多普勒, 实际产生的多普勒是cfg.max_norDoppler-5，累积到雷达图边缘
cfg.max_nor_delay_no_upsample = para.Ncp - 1;
cfg.max_nor_delay = cfg.max_nor_delay_no_upsample * para.L - 1;  % 接收端粗搜索最大归一化时延（过采样后）
cfg.min_nor_delay = 0;                      % 接收端粗搜索最小归一化时延

cfg.min_nor_doppler = -cfg.max_nor_doppler_abs * para.L / 2;
cfg.max_nor_doppler = cfg.max_nor_doppler_abs * para.L / 2 - 1 / para.L_doppler;

cfg.decimals = 5;               % 时延和多普勒保留小数位数
cfg.guard_interval = 5;         % 时延和多普勒的保护间隔
% 归一化时延随机范围
cfg.random_delay_range = [cfg.guard_interval, para.Ncp-cfg.guard_interval];
% 归一化多普勒随机范围
cfg.random_doppler_range = [-cfg.max_nor_doppler+cfg.guard_interval,  cfg.max_nor_doppler-cfg.guard_interval];
cfg.delay_min_gap = 2;             % 多目标时，不同目标归一化时延的最小间隔
cfg.doppler_min_gap = 2;           % 多目标时，不同目标归一化多普勒的最小间隔

cfg.flag_delayMigration = 1;                    % 1: 符号级距离徙动；2: 采样点级距离徙动；3：不徙动
fprintf('\n================ 徙动模式 ===================\n');
fprintf('距离徙动模式 flag_delayMigration       = %s', delayMigrationStr(cfg.flag_delayMigration));
fprintf('\n============================================\n');

% 波形理论性能：慢时间采用含 CP 和前后保护区的实际符号间隔 T_symbol。
cfg.delta_tau = para.T / para.Nc;                                % 时延分辨率 1/B，单位 s；不是上采样后的采样间隔
cfg.delta_fd = 1 / (para.T_symbol * para.Nsym);                  % 多普勒分辨率 1/CPI，单位 Hz；CPI 为相干处理时长
cfg.delta_tau_nor = cfg.delta_tau * para.fs_up;                % 归一化时延分辨率，采用 tau*fs_up 坐标；当前 B=fs 时等于 L
cfg.delta_fd_nor = cfg.delta_fd / para.delta_f;                % 归一化多普勒分辨率，采用 nu=fd/delta_f 坐标；不是搜索网格步长
cfg.delta_range_com = cfg.delta_tau * para.c;                   % 单程传播距离分辨率 c/B，单位 m
cfg.delta_velocity_com = cfg.delta_fd / para.fc * para.c;       % 单程多普勒对应的径向速度分辨率，单位 m/s
cfg.delta_range_rad = cfg.delta_tau * para.c / 2;               % 单站雷达双程测距分辨率 c/(2B)，单位 m
cfg.delta_velocity_rad = cfg.delta_fd / para.fc * para.c / 2;   % 单站雷达双程多普勒对应的径向速度分辨率，单位 m/s
cfg.Max_unamb_range = para.c * para.T / 2;                       % 由子载波间隔决定的单站距离模糊周期 c/(2*delta_f)，单位 m
cfg.Max_ISI_free_range = para.c * para.Tcp / 2;                  % CP 时长能覆盖的最大双程距离，单位 m；假定同步且无其他时延余量占用
cfg.Max_unamb_velocity = para.c / (2 * para.fc * para.T_symbol); % 慢时间采样对应的速度模糊周期，单位 m/s；正负对称区间为 ±该值/2
cfg.Max_ICI_free_velocity = para.B * para.c / (20 * para.Nc * para.fc); % 按 |fd|<=0.1*delta_f 估算的单站低 ICI 速度幅值，单位 m/s；
% 不是严格无 ICI 界限

fprintf('\n========== 雷达波形感知性能参数 ==========\n');
fprintf('归一化时延分辨率 delta_tau_nor         = %.6g (tau*fs_up)\n', cfg.delta_tau_nor);
fprintf('归一化多普勒分辨率 delta_fd_nor        = %.6g (fd/delta_f)\n', cfg.delta_fd_nor);
fprintf('雷达距离分辨率 delta_range_rad          = %.4f m\n', cfg.delta_range_rad);
fprintf('雷达速度分辨率 delta_velocity_rad      = %.4f m/s\n', cfg.delta_velocity_rad);
fprintf('最大不模糊距离 Max_unamb_range         = %.4f m\n', cfg.Max_unamb_range);
fprintf('最大无ISI距离 Max_ISI_free_range       = %.4f m\n', cfg.Max_ISI_free_range);
fprintf('最大不模糊速度 Max_unamb_velocity      = %.4f m/s\n', cfg.Max_unamb_velocity);
fprintf('最大ICI-free速度 Max_ICI_free_velocity = %.4f m/s\n', cfg.Max_ICI_free_velocity);
fprintf('通信单程距离分辨率 delta_range_com      = %.4f m\n', cfg.delta_range_com);
fprintf('通信单程速度分辨率 delta_velocity_com  = %.4f m/s\n', cfg.delta_velocity_com);
fprintf('==========================================\n\n');

%% ===================== 2.参数搜索算法参数设置  =====================
algorithm = struct();

algorithm.coarse_delay_search_step = 1;                           % 粗估计中时延搜索步长（过采样前的步长）
algorithm.coarse_doppler_search_step = 1 / para.L_doppler;        % 粗估计中多普勒搜索步长（过采样前的步长）
% 搜索中的时延步长为algorithm.coarse_delay_search_step，由于过采样对应的网格为algorithm.coarse_delay_search_step_grid
% 为了对应归一化1的时延和多普勒需要的步数
algorithm.coarse_delay_search_step_grid = 1/algorithm.coarse_delay_search_step;
algorithm.coarse_doppler_search_step_grid = 1/algorithm.coarse_doppler_search_step;

% 粗估计中归一化时延和多普勒的搜索范围
algorithm.nor_delay_search_set = cfg.min_nor_delay:algorithm.coarse_delay_search_step:cfg.max_nor_delay;
algorithm.nor_doppler_search_set = cfg.min_nor_doppler:algorithm.coarse_doppler_search_step:cfg.max_nor_doppler;

algorithm.N_delay_peak = 4;             %
algorithm.N_doppler_peak = 4;           %
algorithm.Num_iter_single = 8;          % 单目标最大迭代次数
algorithm.Num_iter_cyclic = 1;          % 多目标最大迭代次数
algorithm.Pfa_target = 1e-3;

%% ===================== 3. 预计算 =====================
% 生成发射端脉冲成形滤波器系数
% 成形滤波器的产生有问题，后续理论模糊函数的产生这里都需要改，未做能量归一化
[g_ps, P_rrc] = f_pulse_shaping_filter(para);

% 模糊函数存储位置及名称
% cachePath = './AmbiguFunc/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v6_20251022.mat';
cachePath = './AmbiguFunc/AmbigFunc_pulsedAFDM_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v1_20261118.mat';
% cachePath = './AmbiguFunc/AmbigFunc_pulsedAFDM_N1024_2Nc1_11_L2_Ldelay2_Ldoppler2_v1_20261118.mat';

ambigData = load(cachePath, 'AmbigFunc_pulsedAFDM_noDelayDop');
AmbigFunc_pulsedAFDM_noDelayDop = ambigData.AmbigFunc_pulsedAFDM_noDelayDop;
AmbigFunc_flip = AmbigFunc_pulsedAFDM_noDelayDop([1, end:-1:2], [1, end:-1:2]);
AmbigFunc_pulsedAFDM_noDelayDop = AmbigFunc_pulsedAFDM_noDelayDop / para.N^2;
AmbigFunc_flip = AmbigFunc_flip / para.N^2;

%% ===================== 4. 感知参数生成 =====================
% 随机生成增益，时延和多普勒
cfg.chan_coef = [1];
cfg.delay_vec = f_generate_random_num_with_gap(cfg.random_delay_range, cfg.num_target, ...
    cfg.decimals, cfg.delay_min_gap);
cfg.doppler_vec = f_generate_random_num_with_gap(cfg.random_doppler_range, cfg.num_target, ...
    cfg.decimals, cfg.doppler_min_gap);

% cfg.doppler_vec = 12.35078;
%% ===================== 5. 发送信号生成 =====================
% 生成该时延和多普勒条件下的目标运动轨迹
[chan_coef, doppler_i_cal, delay_i_cal, tau_vector, true_range, true_velocity_rad, delay_max] = ...
    f_sensing_delay_doppler_trajectory_v5(cfg, para, cfg.delay_vec, cfg.doppler_vec, cfg.chan_coef);

[Tx_signal, Tx_signal_PS_AFDM_eff, ~, ~, ~] = f_AFDM_TxSignal_v2(para, g_ps, P_rrc);

% 每帧只计算一次：由已知发射参考得到主峰功率近似权重，平均值归一化为1。
% 仅用于理论模板及其导数，不给实际接收功率面再次加权。
E_symbol = sum(abs(Tx_signal_PS_AFDM_eff).^2, 1) / (para.N * para.L);
weight_symbol = E_symbol.^2;
weight_symbol = weight_symbol / mean(weight_symbol);
cfg.symbol_power_weight = weight_symbol;
% cfg.symbol_power_weight = ones(size(weight_symbol));
%% ===================== 6. 过信道 =====================
[Rx_echo_signal, ~, power_signal, power_noise] = f_time_vary_channel_v5(...
    Tx_signal(:), para.Nc, para.L, cfg.num_target, chan_coef, tau_vector, ...
    doppler_i_cal, delay_max, para.sinc_filter_len, para.SNR, add_noise_flag);

% 根据噪声设置门限
P_noise_MF = power_noise ./ para.N;
Th_noTarget = 1 / para.Nsym * gaminv(1 - algorithm.Pfa_target, para.Nsym, P_noise_MF);

%% ===================== 7. 感知信号接收 =====================
Rx_PS_remove_GPS = f_AFDM_Rx_preprocess_v1(Rx_echo_signal, cfg, para);

%%%%%%%%%%%%%%%%%%%%%%%%% 收发信号时域实部 %%%%%%%%%%%%%%%%%%%%%%%%
% figure('Name', ['收发信号实部对比']);
% stem(real(Rx_echo_signal(1:2000)), '--r');
% hold on;
% stem(real(Tx_signal(1:2000)), '-bo');
% legend('Rx','Tx');
% % title('收发信号实部对比');
% %%%%%%%%%%%%%%%%%%%%%%%%% 收发信号时域实部 %%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%%%%%%%%%%%%%%%%%%%%%%% 收发信号频谱 %%%%%%%%%%%%%%%%%%%%%%%%
% length_rx = length(Rx_echo_signal);
% f_tx = (-para.len_signal/2:(para.len_signal/2-1))*(para.fs/para.len_signal);
% f_rx = (-length_rx/2:(length_rx/2-1))*(para.fs/length_rx);
%
% Freq_tx_signal = abs(fftshift(fft(Tx_signal)));
% Freq_rx_signal = abs(fftshift(fft(Rx_echo_signal)));
%
% figure('Name', ['收发信号频谱对比']);
% plot(f_rx, 20*log10(Freq_rx_signal/max(Freq_rx_signal)), '-bo');  hold on;
% plot(f_tx, 20*log10(Freq_tx_signal/max(Freq_tx_signal)), '-ro');
% ylim([-150,5]);
% legend('Rx','Tx');
% % title('收发信号频谱对比');
%%%%%%%%%%%%%%%%%%%%%%%%% 收发信号频谱 %%%%%%%%%%%%%%%%%%%%%%%%

%% ===================== 8. 徙动补偿与参数估计算法 =====================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1.RFT算法-start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RFT 频域基线：包络徙动补偿后做慢时间 FFT，读取候选速度对应的频点。
% delay_search_grid = cfg.min_nor_delay:algorithm.coarse_delay_search_step:cfg.max_nor_delay;
% doppler_search_grid = cfg.min_nor_doppler:algorithm.coarse_doppler_search_step:cfg.max_nor_doppler;

% 理论分辨率不等于估计网格步长；每个距离分辨单元取至少 10 个点。
% rft_delay_step = min(algorithm.coarse_delay_search_step_grid, cfg.delta_tau_nor / 10);
% delay_search_grid = cfg.min_nor_delay:rft_delay_step:cfg.max_nor_delay;
% doppler_search_grid = cfg.min_nor_doppler:cfg.delta_fd:cfg.max_nor_doppler;

[delay_search_grid, doppler_search_grid] = f_RFT_searchRange_v1(para.SNR, cfg.delay_vec, ...
    cfg.doppler_vec, para.Nsym, cfg.min_nor_delay, cfg.max_nor_delay, cfg.min_nor_doppler, ...
    cfg.max_nor_doppler, algorithm.coarse_delay_search_step, algorithm.coarse_doppler_search_step);

[~, RFT_noc_paramenter_esti] = f_AFDM_RFT_v5(Rx_PS_remove_GPS, Tx_signal_PS_AFDM_eff, ...
    delay_search_grid, doppler_search_grid, cfg, para);

[~, RFT_no_paramenter_esti_v6, ~, ~] = f_AFDM_RFT_v6(Rx_PS_remove_GPS, ...
    Tx_signal_PS_AFDM_eff, delay_search_grid, doppler_search_grid, cfg, para);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1.RFT算法-end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2.AFNE算法-start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% （1）粗搜索 %%%%%
[ML_value_refine_set_aver, ~, Migration_non, migration_compensate_bin] = ...
    f_AFDM_coarseRM_MF_v12(Rx_PS_remove_GPS, Tx_signal_PS_AFDM_eff, para, algorithm, cfg);

cfg.migration_compensate_bin = migration_compensate_bin;

% 
AFNE_parameters_esti = f_paraEstibyGradDesc_newtonOMP_single_trarget_v9( ...
    ML_value_refine_set_aver, AmbigFunc_flip, algorithm.nor_delay_search_set, algorithm.nor_doppler_search_set, ...
    para.N, para.L, algorithm.coarse_delay_search_step_grid, algorithm.coarse_doppler_search_step_grid, ...
    cfg.max_nor_delay_no_upsample, cfg.max_nor_doppler_abs,algorithm.N_delay_peak, algorithm.N_doppler_peak, ...
    para.fs_up, para.delta_f, algorithm.Num_iter_single, algorithm.Num_iter_cyclic,...
    Th_noTarget, cfg.Max_num_target_est, para, Migration_non, cfg);

fprintf('AFNE归一化时延估计误差 = %.5f \n', cfg.delay_vec - AFNE_parameters_esti(1));
fprintf('AFNE归一化多普勒估计误差 = %.5f \n', cfg.doppler_vec - AFNE_parameters_esti(2));

fprintf('RFT归一化时延估计误差 = %.5f \n', cfg.delay_vec - RFT_no_paramenter_esti_v6(1));
fprintf('RFT归一化多普勒估计误差 = %.5f \n', cfg.doppler_vec - RFT_no_paramenter_esti_v6(2));

fprintf('RFT-old归一化时延估计误差 = %.5f \n', cfg.delay_vec - RFT_noc_paramenter_esti(1));
fprintf('RFT-old归一化多普勒估计误差 = %.5f \n', cfg.doppler_vec - RFT_noc_paramenter_esti(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2.AFNE算法-end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
