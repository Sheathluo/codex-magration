%%%%%%%%%%%%%%%%%%%%%%%%%
% 版本信息：main_20261104_v11_Monte_MagriCom_AFDM_N128
% 函数功能：仿真主函数:1.只跑一次，不仿真多次
% 创建日期：2026年11月4日
% 创建人：LuoTao
%
% 版本信息：main_20261114_v12_single_target_Monte_MagriCom_N128
% 函数功能：仿真主函数:1.改成蒙特卡洛仿真;
%                   2.仅支持单目标;
%                   3.只有AFNE算法;
% 创建日期：2026年11月14日
% 创建人：LuoTao
%
% 版本信息：main_20261114_v12_single_target_Monte_MagriCom_N128
% 函数功能：仿真主函数:1.增加RFT非相干累计算法作为对比
% 创建日期：2026年11月16日
% 创建人：LuoTao
%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

%% ===================== 0. 信号基本参数配置 =====================
para = struct();
% 物理常量
para.c = 3e8;

para.Nc = 1024;                  % 子载波数

% 系统基本参数
para.fc = 11.7e9;
para.lambda = para.c/para.fc;
para.fs_up = 10e6;              % 上采样后采样率
if para.Nc == 1024; para.fs_up = 240e6; end
para.L = 2;                     % 过采样倍数
para.L_delay = para.L;          % 理论模糊函数的时延维过采样倍数
para.L_doppler = para.L;        % 理论模糊函数/接收搜索的多普勒维过采样倍数
para.fs = para.fs_up / para.L;  % 采样率
para.Ts = 1/para.fs_up;         % 采样间隔
para.B = para.fs;               % 带宽
para.sinc_filter_len = 300;     % 插值滤波器阶数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para.Nsym = 200;                % 符号数/帧数
add_noise_flag = true;         % 是否加噪 true
para.SNR = -40:5:10;
% para.SNR = -25;
para.Ntrial = 16;
% 随机种子
% rng('shuffle');
% rng(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 脉冲成形参数
para.beta = 0.35;               % RRC/RC 滚降因子
if para.Nc == 1024; para.beta = 0.25; end % RRC/RC 滚降因子 Nc = 1024
para.N_preSuffix = 20;          % 每个 AFDM 符号前后额外保护点长度
para.N_tapFiler = 20;           % 成形滤波器半边抽头参数

% AFDM 基本参数
% para.Nc = 128;                  % 子载波数
para.c1 = 8 / (2 * para.Nc);    % c1
if para.Nc == 1024; para.c1 = 11 / (2 * para.Nc); end  % c1 Nc = 1024
para.c2 = 0;                    % c2
para.Mod = 16;                  % 16-QAM
para.mu_4 = 1.32;               % 16QAM四阶矩
para.N = para.Nc;
para.M = para.Nc;
para.delta_f = para.B / para.Nc;    % 子载波间隔
para.T = 1 / para.delta_f;          % 有效符号持续时间
para.Ncp = 20;                      % 循环前缀点数(上采样前)
if para.Nc == 1024; para.Ncp = 32; end  % 循环前缀点数(上采样前) Nc = 1024
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
cfg.num_target = 1;                         % 目标数
cfg.Max_num_target_est = cfg.num_target;    % 接收端 OMP 最多检测目标数，不一定等于真实目标数

cfg.max_nor_doppler_abs = 15;               % 接收端粗搜索最大归一化多普勒, 实际产生的多普勒是cfg.max_norDoppler-5，累积到雷达图边缘
cfg.max_nor_delay_no_upsample = para.Ncp - 1;
cfg.max_nor_delay = cfg.max_nor_delay_no_upsample * para.L - 1;  % 接收端粗搜索最大归一化时延（过采样后）
cfg.min_nor_delay = 0;                      % 接收端粗搜索最小归一化时延

cfg.min_nor_doppler = -cfg.max_nor_doppler_abs * para.L / 2;
cfg.max_nor_doppler = cfg.max_nor_doppler_abs * para.L / 2 - 1 / para.L_doppler;

cfg.decimals = 5;               % 时延和多普勒保留小数位数
cfg.guard_interval = 4;         % 时延和多普勒的保护间隔
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
algorithm.Num_iter_cyclic = algorithm.Num_iter_single;          % 多目标最大迭代次数
algorithm.Pfa_target = 1e-3;

%% ===================== 3. 预计算 =====================
% 生成发射端脉冲成形滤波器系数
% 成形滤波器的产生有问题，后续理论模糊函数的产生这里都需要改，未做能量归一化
[g_ps, P_rrc] = f_pulse_shaping_filter(para);

% 模糊函数存储位置及名称
% cachePath = './AmbiguFunc/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v6_20251022.mat';
% cachePath = './AmbiguFunc/AmbigFunc_pulsedAFDM_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v1_20261118.mat';
cachePath = './AmbiguFunc/AmbigFunc_pulsedAFDM_N1024_2Nc1_11_L2_Ldelay2_Ldoppler2_v1_20261118.mat';

ambigData = load(cachePath, 'AmbigFunc_pulsedAFDM_noDelayDop');
AmbigFunc_pulsedAFDM_noDelayDop = ambigData.AmbigFunc_pulsedAFDM_noDelayDop;
AmbigFunc_flip = AmbigFunc_pulsedAFDM_noDelayDop([1, end:-1:2], [1, end:-1:2]);
AmbigFunc_pulsedAFDM_noDelayDop = AmbigFunc_pulsedAFDM_noDelayDop / para.N^2;
AmbigFunc_flip = AmbigFunc_flip / para.N^2;

%% ===================== 4. 蒙特卡洛仿真 =====================
snr_len = numel(para.SNR);
nTrial = para.Ntrial;
% 每次仿真的真实参数
delay_real_matrix = zeros(snr_len, para.Ntrial, cfg.num_target);
doppler_real_matrix = zeros(size(delay_real_matrix));
channel_coef_real_matrix = zeros(size(delay_real_matrix));

% 每次仿真的估计参数
AFNE_delay_esti_matrix = zeros(size(delay_real_matrix));
AFNE_doppler_esti_matrix = zeros(size(delay_real_matrix));
AFNE_channel_coef_esti_matrix = zeros(size(delay_real_matrix));
RFT_delay_esti_matrix = zeros(size(delay_real_matrix));
RFT_doppler_esti_matrix = zeros(size(delay_real_matrix));
RFT_channel_coef_esti_matrix = zeros(size(delay_real_matrix));

RFT_delay_esti_matrix_old = zeros(size(delay_real_matrix));
RFT_doppler_esti_matrix_old = zeros(size(delay_real_matrix));
RFT_channel_coef_esti_matrix_old = zeros(size(delay_real_matrix));

% 每次仿真的参数估计误差
AFNE_error_delay_matrix = zeros(size(delay_real_matrix));
AFNE_error_doppler_matrix = zeros(size(delay_real_matrix));
AFNE_error_channel_coef_matrix = zeros(size(delay_real_matrix));
RFT_error_delay_matrix = zeros(size(delay_real_matrix));
RFT_error_doppler_matrix = zeros(size(delay_real_matrix));
RFT_error_channel_coef_matrix = zeros(size(delay_real_matrix));

RFT_error_delay_matrix_old = zeros(size(delay_real_matrix));
RFT_error_doppler_matrix_old = zeros(size(delay_real_matrix));
RFT_error_channel_coef_matrix_old = zeros(size(delay_real_matrix));

% 每个SNR下估计误差的RMSE
AFNE_delay_esti_RMSE = zeros(1,snr_len);
AFNE_doppler_esti_RMSE = zeros(1,snr_len);
AFNE_channel_coef_esti_RMSE = zeros(1,snr_len);
RFT_delay_esti_RMSE = zeros(1,snr_len);
RFT_doppler_esti_RMSE = zeros(1,snr_len);
RFT_channel_coef_esti_RMSE = zeros(1,snr_len);

RFT_delay_esti_RMSE_old = zeros(1,snr_len);
RFT_doppler_esti_RMSE_old = zeros(1,snr_len);
RFT_channel_coef_esti_RMSE_old = zeros(1,snr_len);

for snr_idx = 1:snr_len
    snr = para.SNR(snr_idx);

    % 每个 worker 独立随机流，保证可复现且互不相关
    % rs = RandStream('Threefry', 'Seed', 1000 + snr_idx);
    tic;                    % 开始计时
    parfor trial_idx = 1:nTrial
    % for trial_idx = 1:nTrial
        cfg_trail = cfg;
        % ===================== 4.1 生成感知参数 =====================
        cfg_trail.chan_coef = [1];
        cfg_trail.delay_vec = f_generate_random_num_with_gap(cfg_trail.random_delay_range, cfg_trail.num_target, ...
            cfg_trail.decimals, cfg_trail.delay_min_gap);
        cfg_trail.doppler_vec = f_generate_random_num_with_gap(cfg_trail.random_doppler_range, cfg_trail.num_target, ...
            cfg_trail.decimals, cfg_trail.doppler_min_gap);
        % 保存每次的时延、多普勒和增益
        delay_real_matrix(snr_idx, trial_idx, :) = cfg_trail.delay_vec;
        doppler_real_matrix(snr_idx, trial_idx, :) = cfg_trail.doppler_vec;
        channel_coef_real_matrix(snr_idx, trial_idx, :) = cfg_trail.chan_coef;
        % ===================== 4.2. 生成发送信号 =====================
        % 生成该时延和多普勒条件下的目标运动轨迹
        [chan_coef, doppler_i_cal, delay_i_cal, tau_vector, true_range, true_velocity_rad, delay_max] = ...
            f_sensing_delay_doppler_trajectory_v5(cfg_trail, para, cfg_trail.delay_vec, ...
            cfg_trail.doppler_vec, cfg_trail.chan_coef);

        [Tx_signal, Tx_signal_PS_AFDM_eff, ~, ~, ~] = f_AFDM_TxSignal_v2(para, g_ps, P_rrc);

        % 每帧只计算一次：由已知发射参考得到主峰功率近似权重，平均值归一化为1。
        % 仅用于理论模板及其导数，不给实际接收功率面再次加权。
        E_symbol = sum(abs(Tx_signal_PS_AFDM_eff).^2, 1) / (para.N * para.L);
        weight_symbol = E_symbol.^2;
        weight_symbol = weight_symbol / mean(weight_symbol);
        cfg_trail.symbol_power_weight = weight_symbol;
        % cfg_trail.symbol_power_weight = ones(size(weight_symbol));

        % ===================== 4.3. 信道 =====================
        [Rx_echo_signal, ~, power_signal, power_noise] = f_time_vary_channel_v5(...
            Tx_signal(:), para.Nc, para.L, cfg_trail.num_target, chan_coef, tau_vector, ...
            doppler_i_cal, delay_max, para.sinc_filter_len, snr, add_noise_flag);

        % 根据噪声设置门限
        P_noise_MF = power_noise ./ para.N;
        Th_noTarget = 1 / para.Nsym * gaminv(1 - algorithm.Pfa_target, para.Nsym, P_noise_MF);

        % ===================== 4.4. 感知参数接收 =====================
        Rx_PS_remove_GPS = f_AFDM_Rx_preprocess_v1(Rx_echo_signal, cfg_trail, para);

        % ===================== 4.5. AFNE徙动补偿与参数估计算法 =====================
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% （1）粗搜索 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ML_value_refine_set_aver, ~, Migration_non, migration_compensate_bin] = ...
            f_AFDM_coarseRM_MF_v12(Rx_PS_remove_GPS, Tx_signal_PS_AFDM_eff, para, algorithm, cfg_trail);

        cfg_trail.migration_compensate_bin = migration_compensate_bin;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% （2）精搜索 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 目前仅单目标，后续再修改为多目标
        AFNE_parameters_esti = f_paraEstibyGradDesc_newtonOMP_single_trarget_v9( ...
            ML_value_refine_set_aver, AmbigFunc_flip, algorithm.nor_delay_search_set, algorithm.nor_doppler_search_set, ...
            para.N, para.L, algorithm.coarse_delay_search_step_grid, algorithm.coarse_doppler_search_step_grid, ...
            cfg_trail.max_nor_delay_no_upsample, cfg_trail.max_nor_doppler_abs,algorithm.N_delay_peak, ...
            algorithm.N_doppler_peak, para.fs_up, para.delta_f, algorithm.Num_iter_single, ...
            algorithm.Num_iter_cyclic, Th_noTarget, cfg_trail.Max_num_target_est, para, Migration_non, cfg_trail);
        % AFNE_parameters_esti = [0,0,0];
        % 估计值
        AFNE_estimation_delay = AFNE_parameters_esti(1);
        AFNE_estimation_doppler = AFNE_parameters_esti(2);
        AFNE_estimation_channel_coef = AFNE_parameters_esti(3);
        % 保存每次仿真的估计值
        AFNE_delay_esti_matrix(snr_idx, trial_idx, :) = AFNE_estimation_delay;
        AFNE_doppler_esti_matrix(snr_idx, trial_idx, :) = AFNE_estimation_doppler;
        AFNE_channel_coef_esti_matrix(snr_idx, trial_idx, :) = AFNE_estimation_channel_coef;
        % 保存每次估计误差
        AFNE_error_delay_matrix(snr_idx, trial_idx, :) = cfg_trail.delay_vec - AFNE_estimation_delay;
        AFNE_error_doppler_matrix(snr_idx, trial_idx, :) = cfg_trail.doppler_vec - AFNE_estimation_doppler;
        AFNE_error_channel_coef_matrix(snr_idx, trial_idx, :) = cfg_trail.chan_coef - AFNE_estimation_channel_coef;

        % ===================== 4.6. RFT算法 =====================
        % 为了减少计算量，根据经验来设置参数搜索范围
        [delay_search_grid, doppler_search_grid] = f_RFT_searchRange_v1(snr, cfg_trail.delay_vec, ...
            cfg_trail.doppler_vec, para.Nsym, cfg_trail.min_nor_delay, cfg_trail.max_nor_delay, cfg_trail.min_nor_doppler, ...
            cfg_trail.max_nor_doppler, algorithm.coarse_delay_search_step, algorithm.coarse_doppler_search_step);

        [~, RFT_noc_paramenter_esti_old] = f_AFDM_RFT_v5(Rx_PS_remove_GPS, Tx_signal_PS_AFDM_eff, ...
            delay_search_grid, doppler_search_grid, cfg_trail, para);
        [~, RFT_noc_paramenter_esti, ~, ~] = f_AFDM_RFT_v6(Rx_PS_remove_GPS, ...
            Tx_signal_PS_AFDM_eff, delay_search_grid, doppler_search_grid, cfg_trail, para);

        % RFT_noc_paramenter_esti = [0,0];
        RFT_estimation_delay = RFT_noc_paramenter_esti(1);
        RFT_estimation_doppler = RFT_noc_paramenter_esti(2);
        RFT_estimation_channel_coef = 1.01;

        RFT_estimation_delay_old = RFT_noc_paramenter_esti_old(1);
        RFT_estimation_doppler_old = RFT_noc_paramenter_esti_old(2);
        RFT_estimation_channel_coef_old = 1.01;

        % 保存每次仿真的估计值
        RFT_delay_esti_matrix(snr_idx, trial_idx, :) = RFT_estimation_delay;
        RFT_doppler_esti_matrix(snr_idx, trial_idx, :) = RFT_estimation_doppler;
        RFT_channel_coef_esti_matrix(snr_idx, trial_idx, :) = RFT_estimation_channel_coef;

        RFT_delay_esti_matrix_old(snr_idx, trial_idx, :) = RFT_estimation_delay_old;
        RFT_doppler_esti_matrix_old(snr_idx, trial_idx, :) = RFT_estimation_doppler_old;
        RFT_channel_coef_esti_matrix_old(snr_idx, trial_idx, :) = RFT_estimation_channel_coef_old;
        % 保存每次估计误差
        RFT_error_delay_matrix(snr_idx, trial_idx, :) = cfg_trail.delay_vec - RFT_estimation_delay;
        RFT_error_doppler_matrix(snr_idx, trial_idx, :) = cfg_trail.doppler_vec - RFT_estimation_doppler;
        RFT_error_channel_coef_matrix(snr_idx, trial_idx, :) = cfg_trail.chan_coef - RFT_estimation_channel_coef;

        RFT_error_delay_matrix_old(snr_idx, trial_idx, :) = cfg_trail.delay_vec - RFT_estimation_delay_old;
        RFT_error_doppler_matrix_old(snr_idx, trial_idx, :) = cfg_trail.doppler_vec - RFT_estimation_doppler_old;
        RFT_error_channel_coef_matrix_old(snr_idx, trial_idx, :) = cfg_trail.chan_coef - RFT_estimation_channel_coef_old;

        if mod(trial_idx, 50) == 0
            fprintf('\nSNR = %d dB的第 %d 轮迭代:\n',snr, trial_idx);
            fprintf('真实归一化时延 = %.5f \n', cfg_trail.delay_vec);
            % fprintf('AFNE估计归一化时延= %.5f \n', AFNE_estimation_delay);
            fprintf('AFNE归一化时延估计误差 = %.5f \n', cfg_trail.delay_vec - AFNE_estimation_delay);
            fprintf('RFT归一化时延估计误差 = %.5f \n', cfg_trail.delay_vec - RFT_estimation_delay);
            fprintf('RFT-old归一化时延估计误差 = %.5f \n', cfg_trail.delay_vec - RFT_estimation_delay_old);

            % fprintf('\n');
            fprintf('真实归一化多普勒 = %.5f \n', cfg_trail.doppler_vec);
            % fprintf('AFNE估计归一化多普勒= %.5f \n', estimation_doppler);
            fprintf('AFNE归一化多普勒估计误差 = %.5f \n', cfg_trail.doppler_vec - AFNE_estimation_doppler);
            fprintf('RFT归一化多普勒估计误差 = %.5f \n', cfg_trail.doppler_vec - RFT_estimation_doppler);
            fprintf('RFT-old归一化多普勒估计误差 = %.5f \n', cfg_trail.doppler_vec - RFT_estimation_doppler_old);

            % fprintf('\n');
            % fprintf('真实幅度= %.5f \n', cfg_trail.chan_coef);
            % fprintf('估计幅度= %.5f \n', estimation_channel_coef);
            % fprintf('幅度估计误差 = %.5f \n', cfg_trail.chan_coef - estimation_channel_coef);
        end
    end

    elapsedTime = toc;      % 结束计时，返回秒数
    % fprintf('\n');
    
    % 计划运行的蒙特卡洛仿真次数
    % nTrial_ = 1500;
    % Estimated_time = elapsedTime * (nTrial_ / nTrial) * snr_len / 3600;
    % fprintf('仿真%d次预计耗时: %.2f 小时\n', nTrial_, Estimated_time);

    AFNE_delay_esti_RMSE = sqrt(mean(AFNE_error_delay_matrix.^2, 2));
    AFNE_doppler_esti_RMSE = sqrt(mean(AFNE_error_doppler_matrix.^2, 2));
    AFNE_channel_coef_esti_RMSE = sqrt(mean(AFNE_error_channel_coef_matrix.^2, 2));

    RFT_delay_esti_RMSE = sqrt(mean(RFT_error_delay_matrix.^2, 2));
    RFT_doppler_esti_RMSE = sqrt(mean(RFT_error_doppler_matrix.^2, 2));
    RFT_channel_coef_esti_RMSE = sqrt(mean(RFT_error_channel_coef_matrix.^2, 2));

    RFT_delay_esti_RMSE_old = sqrt(mean(RFT_error_delay_matrix_old.^2, 2));
    RFT_doppler_esti_RMSE_old = sqrt(mean(RFT_error_doppler_matrix_old.^2, 2));
    RFT_channel_coef_esti_RMSE_old = sqrt(mean(RFT_error_channel_coef_matrix_old.^2, 2));

    fprintf('\nSNR=%d dB仿真完成，耗时: %.4f 秒', snr, elapsedTime);
    fprintf('\n AFNE delay RMSE = %.5f', AFNE_delay_esti_RMSE(snr_idx));
    fprintf('\n AFNE doppler RMSE = %.5f', AFNE_doppler_esti_RMSE(snr_idx));

    fprintf('\n RFT delay RMSE = %.5f', RFT_delay_esti_RMSE(snr_idx));
    fprintf('\n RFT doppler RMSE = %.5f', RFT_doppler_esti_RMSE(snr_idx));

    fprintf('\n old RFT delay RMSE = %.5f', RFT_delay_esti_RMSE_old(snr_idx));
    fprintf('\n old RFT doppler RMSE = %.5f \n', RFT_doppler_esti_RMSE_old(snr_idx));
end


figure;
semilogy(para.SNR, AFNE_delay_esti_RMSE, '-o', 'LineWidth', 1.5);
hold on;
semilogy(para.SNR, RFT_delay_esti_RMSE, 'r+-', 'LineWidth', 1.5);
hold on;
semilogy(para.SNR, RFT_delay_esti_RMSE_old, 'g+-', 'LineWidth', 1.5);
grid on;
legend('AFNE','RFT','RFT-old');
xlabel('SNR (dB)');
ylabel('Normalized Delay RMSE');
title('AFNE Delay RMSE vs SNR');

figure;
semilogy(para.SNR, AFNE_doppler_esti_RMSE, '-o', 'LineWidth', 1.5);
hold on;
semilogy(para.SNR, RFT_doppler_esti_RMSE, 'r+-', 'LineWidth', 1.5);
hold on;
semilogy(para.SNR, RFT_doppler_esti_RMSE_old, 'g+-', 'LineWidth', 1.5);
grid on;
legend('AFNE','RFT','RFT-old');
xlabel('SNR (dB)');
ylabel('Normalized Doppler RMSE');
title('AFNE Doppler RMSE vs SNR');


% ===== 自动按日期保存仿真结果 =====
resultRoot = './result';
dateFolder = string(datetime('today', 'Format', 'yyyy-MM-dd'));
timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
% 完整路径
saveFolder = fullfile(resultRoot, dateFolder);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
% cfg.monteCarloResultPath = fullfile(saveFolder, sprintf('AFDM_MC_SNR_RMSE_N128_%s.mat', timestamp));
% ===== 保存为 .mat 文件 =====
matFileName = sprintf('RMSE_result_%s.mat', timestamp);
matFilePath = fullfile(saveFolder, matFileName);

save(matFilePath, ...
    'para', 'cfg', 'algorithm', 'delay_real_matrix', 'doppler_real_matrix', 'channel_coef_real_matrix', ...
    'AFNE_delay_esti_matrix', 'AFNE_doppler_esti_matrix','AFNE_channel_coef_esti_matrix', ...
    'AFNE_delay_esti_RMSE', 'AFNE_doppler_esti_RMSE','AFNE_channel_coef_esti_RMSE', ...
    'RFT_delay_esti_matrix', 'RFT_doppler_esti_matrix', 'RFT_channel_coef_esti_matrix', ...
    'RFT_delay_esti_RMSE', 'RFT_doppler_esti_RMSE','RFT_channel_coef_esti_RMSE',...
    'RFT_delay_esti_matrix_old', 'RFT_doppler_esti_matrix_old', 'RFT_channel_coef_esti_matrix_old', ...
    'RFT_delay_esti_RMSE_old', 'RFT_doppler_esti_RMSE_old','RFT_channel_coef_esti_RMSE_old');


% matFilePath = '/Users/luotao/Documents/Codex/codex-magration/result/2026-11-19/RMSE_result_20261119_161455.mat';
% matFilePath = 'E:\Codex\codex-magration\result\2026-11-20\RMSE_result_20261120_105126.mat';
% load(matFilePath);   % 变量会直接恢复到当前工作区
% figure;
% semilogy(para.SNR, AFNE_delay_esti_RMSE, '-o', 'LineWidth', 1.5);
% hold on;
% semilogy(para.SNR, RFT_delay_esti_RMSE, 'r+-', 'LineWidth', 1.5);
% hold on;
% semilogy(para.SNR, RFT_delay_esti_RMSE_old, 'g+-', 'LineWidth', 1.5);
% grid on;
% legend('AFNE','RFT','RFT-old');
% xlabel('SNR (dB)');
% ylabel('Normalized Delay RMSE');
% title('AFNE Delay RMSE vs SNR');
% 
% figure;
% semilogy(para.SNR, AFNE_doppler_esti_RMSE, '-o', 'LineWidth', 1.5);
% hold on;
% semilogy(para.SNR, RFT_doppler_esti_RMSE, 'r+-', 'LineWidth', 1.5);
% hold on;
% semilogy(para.SNR, RFT_doppler_esti_RMSE_old, 'g+-', 'LineWidth', 1.5);
% grid on;
% legend('AFNE','RFT','RFT-old');
% xlabel('SNR (dB)');
% ylabel('Normalized Doppler RMSE');
% title('AFNE Doppler RMSE vs SNR');
