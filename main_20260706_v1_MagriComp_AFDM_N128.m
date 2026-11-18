%% AFDM 蒙特卡洛 SNR 扫描：发射 - 时变感知信道 - 接收估计
% 主脚本保留实验配置、Monte Carlo 主循环、统计、画图和保存。
% 发射、接收预处理、匹配滤波等链路模块由 f_* 函数承接。
% Magrition Compensate
clear;
close all;
tic;

%% ================== 显示辅助函数 ==================
function str = delayMigrationStr(x)
    switch x
        case 1
            str = '符号级徙动';
        case 2
            str = '采样点级徙动';
        case 3
            str = '不徙动';
        otherwise
            str = '未知模式';
    end
end

onoff = @(x) (x ~= 0) * '开' + (x == 0) * '关';

%% ===================== 0. 参数配置 =====================
cfg = struct();

% 蒙特卡洛控制
cfg.enablePlots = true;                        % 保留兼容字段；parfor 内不画图，诊断图由 enableDiagnosticPlots 控制
cfg.addNoiseFlag = false;                        % 是否加噪
% cfg.SNR_dB_set = -10:5:30;                      % 蒙特卡洛扫描的 SNR 点
cfg.SNR_dB_set = 30;                           % 蒙特卡洛扫描的 SNR 点
cfg.numMonteCarlo = 1;                        % 每个SNR点的蒙特卡洛次数
cfg.enableSummaryPlots = true;                  % 是否绘制 RMSE-vs-SNR 汇总曲线
cfg.saveMonteCarloResult = true;                % 是否保存蒙特卡洛原始结果和统计量
cfg.enableDiagnosticPlots = cfg.enablePlots;    % 是否绘制指定 trial 的匹配滤波/残差诊断图
cfg.diagnosticSNRIndex = 1;                     % 诊断图对应的 SNR 索引
cfg.diagnosticTrialIndex = 1;                   % 诊断图对应的 trial 索引

fprintf('\n================ 仿真控制 ===================\n');
fprintf('画图标志位 enablePlots                    = %s\n', onoff(cfg.enablePlots));
fprintf('是否加噪 addNoiseFlag                   = %s\n', onoff(cfg.addNoiseFlag));
fprintf('SNR 扫描范围 SNR_dB_set = %s\n', mat2str(cfg.SNR_dB_set));
fprintf('蒙特卡洛次数 numMonteCarlo             = %d\n', cfg.numMonteCarlo);
fprintf('汇总曲线开关 enableSummaryPlots        = %s\n', onoff(cfg.enableSummaryPlots));
fprintf('保存结果开关 saveMonteCarloResult      = %s\n', onoff(cfg.saveMonteCarloResult));
fprintf('诊断图开关 enableDiagnosticPlots       = %s\n', onoff(cfg.enableDiagnosticPlots));
if cfg.enablePlots
    fprintf('诊断图 SNR 索引 diagnosticSNRIndex      = %d\n', cfg.diagnosticSNRIndex);
    fprintf('诊断图 Trial 索引 diagnosticTrialIndex  = %d', cfg.diagnosticTrialIndex);
end
fprintf('\n============================================\n');

% 接收截取与搜索
cfg.index_start = 0;           % 接收信号起始截取点
cfg.rx_L_delay = 1;            % 接收端时延搜索步进为 1/rx_L_delay
cfg.max_norDoppler = 20;       % 接收端粗搜索最大归一化多普勒

% 缓存与输出
cfg.ambigCachePath = './AmbiguFunc/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v6_20251022.mat';
% cfg.monteCarloResultPath = './result/AFDM_MC_SNR_RMSE_N128_20260616.mat';
% ===== 自动按日期保存仿真结果 =====
resultRoot = './result';
dateFolder = string(datetime('today', 'Format', 'yyyy-MM-dd'));
timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
% 完整路径
saveFolder = fullfile(resultRoot, dateFolder);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
cfg.monteCarloResultPath = fullfile(saveFolder, ...
    sprintf('AFDM_MC_SNR_RMSE_N128_%s.mat', timestamp));

% rng('shuffle');                % 蒙特卡洛仿真不固定随机种子
rng(1);
para = struct();

% 脉冲成形参数
para.L = 2;                    % 过采样倍数
para.beta = 0.35;              % RRC/RC 滚降因子
para.L_delay = para.L;         % 理论模糊函数的时延维过采样倍数
para.L_Doppler = para.L;       % 理论模糊函数/接收搜索的多普勒维过采样倍数
para.N_preSuffix = 20;         % 每个 AFDM 符号前后额外保护点长度
para.N_tapFiler = 20;          % 成形滤波器半边抽头参数
para.fs_up = 80e6 / 8;         % 上采样后采样率
para.fs = para.fs_up / para.L;

% AFDM 基本参数
para.c = 3e8;
para.fc = 24e9;
para.Nc = 128;                 % 子载波数
para.N = para.Nc;
para.M = para.Nc;
para.Nsym = 300;               % 符号数/帧数
para.Mod = 16;                 % 16-QAM
para.B = para.fs;
para.delta_f = para.B / para.Nc;
para.T = 1 / para.delta_f;
para.Ncp = 20;                 % 循环前缀点数
para.Tcp = para.Ncp / para.B;
para.T_AFDM = para.T + para.Tcp;
para.c1 = 8 / (2 * para.Nc);   % AFDM chirp 参数
para.c2 = 0;
para.mu_4 = 1.32;

fprintf('\n=============== 信号基本物理量 ==============\n');
fprintf('信号带宽 B         = %.1f MHz\n', para.B/1e6);
fprintf('采样率  fs         = %.1f Msps\n', para.fs/1e6);
fprintf('累积时间  fs         = %.4f s\n', para.T_AFDM*para.Nsym);
fprintf('\n============================================\n');


% 感知信道目标参数
cfg.Max_num_target_est = 2;    % 接收端 OMP 最多检测目标数，不一定等于真实目标数
cfg.chan_coef = [1];             % 目标/路径复增益

% 每个 Monte Carlo trial 随机生成目标归一化时延/多普勒。
% 这里配置的是随机范围和目标间最小间隔；真实 delay/Doppler 会在 trial 内生成。
cfg.randomParamDecimals = 4;     % 时延和多普勒保留小数位数
cfg.randomDelayRange = [0.5, 25];      % 归一化时延随机范围
cfg.randomDopplerRange = [-cfg.max_norDoppler, cfg.max_norDoppler];    % 归一化多普勒随机范围
cfg.randomDelayMinGap = 2;             % 多目标时，不同目标归一化时延的最小间隔
cfg.randomDopplerMinGap = 2;           % 多目标时，不同目标归一化多普勒的最小间隔

cfg.sincFilterLen = 300;       % 时变分数时延 sinc 插值滤波器长度
cfg.flag_delayMigration = 1;   % 1: 符号级距离徙动；2: 采样点级距离徙动；3：不徙动
fprintf('\n================ 徙动模式 ===================\n');
fprintf('距离徙动模式 flag_delayMigration       = %s', delayMigrationStr(cfg.flag_delayMigration));
fprintf('\n============================================\n');

cfg.chan_coef = cfg.chan_coef(:);
% 目标个数由复增益个数决定；每个 trial 会为这些目标随机分配 delay/Doppler。
cfg.Num_target = length(cfg.chan_coef);
if cfg.Num_target < 1
    error('cfg.chan_coef must contain at least one target/path gain.');
end

% 雷达性能
para.delta_tau = para.T / para.Nc;
para.delta_fd = 1 / (para.T_AFDM * para.Nsym);
para.delta_range_com = para.delta_tau * para.c;
para.delta_velocity_com = para.delta_fd / para.fc * para.c;
para.delta_range_rad = para.delta_tau * para.c / 2;
para.delta_velocity_rad = para.delta_fd / para.fc * para.c / 2;
para.Max_unamb_range = para.c * para.T / 2;
para.Max_ISI_free_range = para.c * para.Tcp / 2;
para.Max_unamb_velocity = para.c / (2 * para.fc * para.T_AFDM);
para.Max_ICI_free_velocity = para.B * para.c / (20 * para.Nc * para.fc);

fprintf('\n========== 雷达波形感知性能参数 ==========\n');
fprintf('雷达距离分辨率 delta_range_rad          = %.4f m\n', para.delta_range_rad);
fprintf('雷达速度分辨率 delta_velocity_rad      = %.4f m/s\n', para.delta_velocity_rad);
fprintf('最大不模糊距离 Max_unamb_range         = %.4f m\n', para.Max_unamb_range);
fprintf('最大无ISI距离 Max_ISI_free_range       = %.4f m\n', para.Max_ISI_free_range);
fprintf('最大不模糊速度 Max_unamb_velocity      = %.4f m/s\n', para.Max_unamb_velocity);
fprintf('最大ICI-free速度 Max_ICI_free_velocity = %.4f m/s\n', para.Max_ICI_free_velocity);
fprintf('通信单程距离分辨率 delta_range_com      = %.4f m\n', para.delta_range_com);
fprintf('通信单程速度分辨率 delta_velocity_com  = %.4f m/s\n', para.delta_velocity_com);
fprintf('==========================================\n\n');

%% ===================== 1. 固定预计算 =====================
% 一个脉冲成形 AFDM 符号包含有效子载波、CP 和前后保护采样。
samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;
Len_sig = samples_per_symbol * para.L * para.Nsym + para.N_preSuffix * para.L;
mod_para = struct('Mod', para.Mod, 'Tb', 1/20/1e6, 'fs', 20e6, 'Rb', 20e6);

% 生成发射端脉冲成形滤波器系数，以及无 CP 有效 AFDM 符号的匹配滤波参考矩阵。
[g_ps, P_rrc] = f_AFDM_pulseShaping_matrix_v1(para);

factor_delay_refine = 1 / cfg.rx_L_delay;
factor_Doppler_refine = 1 / para.L_Doppler;

max_norDelay = para.Ncp - 1;
Min_delay = 0;
Max_delay = max_norDelay * para.L - 1;
Min_nor_Doppler = -cfg.max_norDoppler * para.L / 2;
Max_nor_Doppler = cfg.max_norDoppler * para.L / 2 - 1 / para.L_Doppler;

Nor_delay_Search_set = Min_delay:factor_delay_refine:Max_delay;
Nor_Doppler_Search_set = Min_nor_Doppler:factor_Doppler_refine:Max_nor_Doppler;

% 理论模糊函数用于精估计阶段的目标回波重构/扣除；函数内部会优先读取缓存。
[AmbigFunc_pulsedAFDM_noDelayDop, AmbigFunc_flip] = f_getAmbigFunc_pulseShaped_AFDM_v1(cfg.ambigCachePath);
AmbigFunc_pulsedAFDM_noDelayDop = AmbigFunc_pulsedAFDM_noDelayDop / para.N^2;
AmbigFunc_flip = AmbigFunc_flip / para.N^2;

% OMP/梯度精估计相关参数：峰值搜索窗口、单目标迭代和循环迭代次数、检测虚警率。
N_delay_peak = 4;
N_Doppler_peak = 4;
Num_iter_single = 1;
Num_iter_cyclic = 30;
Pfa_target = 1e-3;

%% ===================== 2. 蒙特卡洛 状态初始化 =====================
SNR_run_set = cfg.SNR_dB_set(:).';
num_trial_run = cfg.numMonteCarlo;
num_snr = length(SNR_run_set);
if cfg.enableDiagnosticPlots
    if cfg.diagnosticSNRIndex < 1 || cfg.diagnosticSNRIndex > num_snr || ...
            cfg.diagnosticTrialIndex < 1 || cfg.diagnosticTrialIndex > num_trial_run
        error('Diagnostic indices must satisfy 1 <= cfg.diagnosticSNRIndex <= %d and 1 <= cfg.diagnosticTrialIndex <= %d.', ...
            num_snr, num_trial_run);
    end
end

% 预分配 trial 级结果矩阵；NaN 表示该 trial 没有产生有效估计。
nor_delay_esti_set = nan(num_snr, num_trial_run);
nor_Doppler_esti_set = nan(num_snr, num_trial_run);
P_sig_strong_esti_set = nan(num_snr, num_trial_run);
range_esti_rad_set = nan(num_snr, num_trial_run);
velocity_esti_rad_set = nan(num_snr, num_trial_run);
nor_delay_error = nan(num_snr, num_trial_run);
nor_Doppler_error = nan(num_snr, num_trial_run);
range_error_m = nan(num_snr, num_trial_run);
velocity_error_mps = nan(num_snr, num_trial_run);
detect_failed = false(num_snr, num_trial_run);
trial_message = cell(num_snr, num_trial_run);

true_nor_delay_set = nan(num_snr, num_trial_run, cfg.Num_target);
true_nor_Doppler_set = nan(num_snr, num_trial_run, cfg.Num_target);
true_range_rad_set = nan(num_snr, num_trial_run, cfg.Num_target);
true_velocity_rad_set = nan(num_snr, num_trial_run, cfg.Num_target);

% *_ref_* 保存当前统计口径使用的参考目标：真实最强目标。
% 多目标时暂不做目标配对，仍用估计端第一列与该参考目标计算 RMSE。
ref_target_index_set = nan(num_snr, num_trial_run);
true_nor_delay_ref_set = nan(num_snr, num_trial_run);
true_nor_Doppler_ref_set = nan(num_snr, num_trial_run);
true_range_ref_rad_set = nan(num_snr, num_trial_run);
true_velocity_ref_rad_set = nan(num_snr, num_trial_run);

nor_delay_esti = NaN;
nor_Doppler_esti = NaN;
P_sig_strong_esti = NaN;
range_esti_rad = NaN;
velocity_esti_rad = NaN;
true_nor_delay = NaN;
true_nor_Doppler = NaN;
true_range_ref_rad = NaN;
true_velocity_ref_rad = NaN;
idx_ref_target = NaN;

%% ===================== 3. SNR 扫描 + 蒙特卡洛 接收估计 =====================
fprintf('Running AFDM 蒙特卡洛 SNR sweep: %d SNR points x %d trials.\n', num_snr, num_trial_run);
progress_queue = parallel.pool.DataQueue;
afterEach(progress_queue, @(msg) fprintf('%s\n', msg));
for idx_snr = 1:num_snr
    SNR_dB_curr = SNR_run_set(idx_snr);
    fprintf('Starting SNR %+g dB: %d trials\n', SNR_dB_curr, num_trial_run);

    nor_delay_esti_row = nan(1, num_trial_run);
    nor_Doppler_esti_row = nan(1, num_trial_run);
    P_sig_strong_esti_row = nan(1, num_trial_run);
    range_esti_rad_row = nan(1, num_trial_run);
    velocity_esti_rad_row = nan(1, num_trial_run);
    nor_delay_error_row = nan(1, num_trial_run);
    nor_Doppler_error_row = nan(1, num_trial_run);
    range_error_m_row = nan(1, num_trial_run);
    velocity_error_mps_row = nan(1, num_trial_run);
    detect_failed_row = false(1, num_trial_run);
    trial_message_row = cell(1, num_trial_run);
    diagnostic_ML_value_row = cell(1, num_trial_run);
    diagnostic_parameters_row = cell(1, num_trial_run);
    diagnostic_true_parameters_row = cell(1, num_trial_run);

    % 当前 SNR 下的 trial 局部缓存，parfor/for 结束后再统一写回全局矩阵。
    true_nor_delay_row = nan(num_trial_run, cfg.Num_target);
    true_nor_Doppler_row = nan(num_trial_run, cfg.Num_target);
    true_range_rad_row = nan(num_trial_run, cfg.Num_target);
    true_velocity_rad_row = nan(num_trial_run, cfg.Num_target);
    ref_target_index_row = nan(1, num_trial_run);
    true_nor_delay_ref_row = nan(1, num_trial_run);
    true_nor_Doppler_ref_row = nan(1, num_trial_run);
    true_range_ref_rad_row = nan(1, num_trial_run);
    true_velocity_ref_rad_row = nan(1, num_trial_run);

    % parfor idx_trial = 1:num_trial_run
    for idx_trial = 1:num_trial_run
        trial_timer = tic;
        try
            % 每个 trial 独立随机目标参数，并据此生成时变时延轨迹和真实误差基准。
            cfg_trial = cfg;
            % cfg_trial.delay_i_cal = f_generate_random_num_with_gap( ...
            %     cfg.randomDelayRange, cfg.Num_target, cfg.randomParamDecimals, cfg.randomDelayMinGap);
            cfg_trial.delay_i_cal = 10;
            % cfg_trial.doppler_i_cal = f_generate_random_num_with_gap( ...
            %     cfg.randomDopplerRange, cfg.Num_target, cfg.randomParamDecimals, cfg.randomDopplerMinGap);
            cfg_trial.doppler_i_cal = -14.5;
            fprintf('Real delay: %.4f \n', cfg_trial.delay_i_cal);
            fprintf('Real doppler: %.4f \n', cfg_trial.doppler_i_cal);
            %%%%%%%%%%%%%%%%%%%% 近似误差计算 %%%%%%%%%%%%%%%%%%%%%%
            true_velocity = cfg_trial.doppler_i_cal * para.delta_f / para.fc * para.c / 2;
            bias = true_velocity/para.fs*(para.Nc*para.L+1)/para.c;         % 单位:s
            bias_sample = bias*para.fs;
            fprintf('bias_sample: %.4f \n', bias_sample);
            %%%%%%%%%%%%%%%%%%%% 近似误差计算 %%%%%%%%%%%%%%%%%%%%%%
            cfg_trial.v_i_cal = cfg_trial.doppler_i_cal;

            % 由本 trial 的随机 delay/Doppler 生成时变时延轨迹，同时得到真实距离/速度。
            [chan_coef, v_i_cal, delay_i_cal, tau_vector, delay_max, true_range_rad, true_velocity_rad] = ...
                f_sensing_delayDoppler_trajectory_v1(cfg_trial, para, samples_per_symbol);
            real_magrition = tau_vector/para.delta_tau;
            % 多目标场景下只统计真实最强目标；估计端仍使用第一列检测结果作为对比对象。
            [~, idx_ref_target_curr] = max(abs(chan_coef));
            true_nor_delay_curr_all = delay_i_cal(:).';
            true_nor_Doppler_curr_all = v_i_cal(:).';
            true_range_rad_curr_all = true_range_rad(:).';
            true_velocity_rad_curr_all = true_velocity_rad(:).';
            true_nor_delay_ref_curr = delay_i_cal(idx_ref_target_curr);
            true_nor_Doppler_ref_curr = v_i_cal(idx_ref_target_curr);
            true_range_ref_rad_curr = true_range_rad(idx_ref_target_curr);
            true_velocity_ref_rad_curr = true_velocity_rad(idx_ref_target_curr);

            % 保存所有真实目标参数，方便后续复盘每个 trial 的随机场景。
            true_nor_delay_row(idx_trial, :) = true_nor_delay_curr_all;
            true_nor_Doppler_row(idx_trial, :) = true_nor_Doppler_curr_all;
            true_range_rad_row(idx_trial, :) = true_range_rad_curr_all;
            true_velocity_rad_row(idx_trial, :) = true_velocity_rad_curr_all;
            ref_target_index_row(idx_trial) = idx_ref_target_curr;
            true_nor_delay_ref_row(idx_trial) = true_nor_delay_ref_curr;
            true_nor_Doppler_ref_row(idx_trial) = true_nor_Doppler_ref_curr;
            true_range_ref_rad_row(idx_trial) = true_range_ref_rad_curr;
            true_velocity_ref_rad_row(idx_trial) = true_velocity_ref_rad_curr;

            % 每个 trial 重新生成随机 QAM 数据、AFDM 发射波形和匹配滤波参考信号。
            [tx_signal, Tx_signal_pulsedAFDM_eff, ~, ~, ~] = ...
                f_AFDM_TxSignal_generator_v1(para, mod_para, g_ps, P_rrc);

            % 时变信道函数负责非整数时延、多普勒、距离徙动，并按当前 SNR 在内部加噪。
            % v4 的多普勒输入是归一化多普勒 nu = f_d / delta_f，因此这里传 v_i_cal 而不是物理速度。
            [Rx_signal_AFDM_tmp, ~, ~, P_noise_cal] = f_time_vary_channel_v4( ...
                tx_signal(:), para.Nc, para.L, cfg.Num_target, chan_coef, tau_vector, ...
                v_i_cal, delay_max, cfg.sincFilterLen, SNR_dB_curr, cfg.addNoiseFlag);

            % 接收预处理：截取有效接收段，去掉脉冲成形过渡带，并提取去保护区后的 AFDM 符号。
            X_pulseShaped_remove_preSuffix = f_AFDM_Rx_preprocess_v1(Rx_signal_AFDM_tmp, cfg, para, Len_sig);

            % 粗距离徙动补偿后，对每个 AFDM 符号做时域匹配滤波并跨符号平均。
            % [ML_value_refine_set_aver, est_dop] = f_AFDM_coarseRM_MF_v3(X_pulseShaped_remove_preSuffix, ...
            %     Tx_signal_pulsedAFDM_eff, para, Nor_Doppler_Search_set, factor_delay_refine, ...
            %     factor_Doppler_refine, Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler);

            [ML_value_refine_set_aver, est_dop] = f_AFDM_coarseRM_MF_v6(X_pulseShaped_remove_preSuffix, ...
                Tx_signal_pulsedAFDM_eff, para, Nor_Doppler_Search_set, factor_delay_refine, ...
                factor_Doppler_refine, Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler);
            fprintf('est_dop: %.4f \n', est_dop);
            % 用 v4 返回的噪声功率设置检测门限，使门限和实际加噪功率口径一致。
            P_noise_MF = P_noise_cal ./ para.N;
            Th_noTarget = 1 / para.Nsym * gaminv(1 - Pfa_target, para.Nsym, P_noise_MF);

            % 基于匹配滤波面和理论模糊函数做 OMP/梯度精估计，输出归一化时延、多普勒和强目标功率。
            set_nor_parameters_esti = f_paraEstibyGradDesc_newtonOMP_multiTarget_v2( ...
                ML_value_refine_set_aver, AmbigFunc_flip, Nor_delay_Search_set, Nor_Doppler_Search_set, ...
                para.N, para.L, cfg.rx_L_delay, para.L_Doppler, max_norDelay, cfg.max_norDoppler, ...
                N_delay_peak, N_Doppler_peak, para.fs_up, para.delta_f, Num_iter_single, Num_iter_cyclic,...
                Th_noTarget, cfg.Max_num_target_est);

            if isempty(set_nor_parameters_esti)
                detect_failed_row(idx_trial) = true;
                trial_message_row{idx_trial} = 'No target was detected.';
                send(progress_queue, sprintf('[no target] SNR %+g dB | trial %d/%d | elapsed %.1f s', ...
                    SNR_dB_curr, idx_trial, num_trial_run, toc(trial_timer)));
            else
                nor_delay_esti_curr = round(set_nor_parameters_esti(1, 1), cfg.randomParamDecimals);
                nor_Doppler_esti_curr = round(set_nor_parameters_esti(2, 1), cfg.randomParamDecimals);

                P_sig_strong_esti_curr = round(set_nor_parameters_esti(3, 1), cfg.randomParamDecimals);
                range_esti_rad_curr = nor_delay_esti_curr / para.fs_up * para.c / 2;
                velocity_esti_rad_curr = nor_Doppler_esti_curr * para.delta_f / para.fc * para.c / 2;

                % 记录当前 trial 的原始估计值和相对真实目标参数的误差
                nor_delay_esti_row(idx_trial) = nor_delay_esti_curr;
                nor_Doppler_esti_row(idx_trial) = nor_Doppler_esti_curr;
                P_sig_strong_esti_row(idx_trial) = P_sig_strong_esti_curr;
                range_esti_rad_row(idx_trial) = range_esti_rad_curr;
                velocity_esti_rad_row(idx_trial) = velocity_esti_rad_curr;
                nor_delay_error_row(idx_trial) = nor_delay_esti_curr - true_nor_delay_ref_curr;
                nor_Doppler_error_row(idx_trial) = nor_Doppler_esti_curr - true_nor_Doppler_ref_curr;
                range_error_m_row(idx_trial) = range_esti_rad_curr - true_range_ref_rad_curr;
                velocity_error_mps_row(idx_trial) = velocity_esti_rad_curr - true_velocity_ref_rad_curr;
                send(progress_queue, sprintf(['[done] SNR %+g dB | trial %d/%d | true delay %.4f | true Doppler %.4f ' ...
                    '| est delay %.4f | est Doppler %.4f | range %.3f m | velocity %.3f m/s | elapsed %.1f s'], ...
                    SNR_dB_curr, idx_trial, num_trial_run, true_nor_delay_ref_curr, true_nor_Doppler_ref_curr, ...
                    nor_delay_esti_curr, nor_Doppler_esti_curr, range_esti_rad_curr, velocity_esti_rad_curr, toc(trial_timer)));

                if cfg.enableDiagnosticPlots && idx_snr == cfg.diagnosticSNRIndex && idx_trial == cfg.diagnosticTrialIndex
                    diagnostic_ML_value_row{idx_trial} = ML_value_refine_set_aver;
                    diagnostic_parameters_row{idx_trial} = [nor_delay_esti_curr; nor_Doppler_esti_curr; P_sig_strong_esti_curr];
                    diagnostic_true_parameters_row{idx_trial} = [true_nor_delay_ref_curr; true_nor_Doppler_ref_curr; ...
                        true_range_ref_rad_curr; true_velocity_ref_rad_curr; idx_ref_target_curr];
                end
            end

        catch ME
            detect_failed_row(idx_trial) = true;
            trial_message_row{idx_trial} = ME.message;
            send(progress_queue, sprintf('[failed] SNR %+g dB | trial %d/%d | %s | elapsed %.1f s', ...
                SNR_dB_curr, idx_trial, num_trial_run, ME.message, toc(trial_timer)));
        end
    end

    nor_delay_esti_set(idx_snr, :) = nor_delay_esti_row;
    nor_Doppler_esti_set(idx_snr, :) = nor_Doppler_esti_row;
    P_sig_strong_esti_set(idx_snr, :) = P_sig_strong_esti_row;
    range_esti_rad_set(idx_snr, :) = range_esti_rad_row;
    velocity_esti_rad_set(idx_snr, :) = velocity_esti_rad_row;
    nor_delay_error(idx_snr, :) = nor_delay_error_row;
    nor_Doppler_error(idx_snr, :) = nor_Doppler_error_row;
    range_error_m(idx_snr, :) = range_error_m_row;
    velocity_error_mps(idx_snr, :) = velocity_error_mps_row;
    detect_failed(idx_snr, :) = detect_failed_row;
    trial_message(idx_snr, :) = trial_message_row;

    % 写回本 SNR 下所有 trial 的真实随机目标参数。
    true_nor_delay_set(idx_snr, :, :) = reshape(true_nor_delay_row, [1, num_trial_run, cfg.Num_target]);
    true_nor_Doppler_set(idx_snr, :, :) = reshape(true_nor_Doppler_row, [1, num_trial_run, cfg.Num_target]);
    true_range_rad_set(idx_snr, :, :) = reshape(true_range_rad_row, [1, num_trial_run, cfg.Num_target]);
    true_velocity_rad_set(idx_snr, :, :) = reshape(true_velocity_rad_row, [1, num_trial_run, cfg.Num_target]);
    ref_target_index_set(idx_snr, :) = ref_target_index_row;
    true_nor_delay_ref_set(idx_snr, :) = true_nor_delay_ref_row;
    true_nor_Doppler_ref_set(idx_snr, :) = true_nor_Doppler_ref_row;
    true_range_ref_rad_set(idx_snr, :) = true_range_ref_rad_row;
    true_velocity_ref_rad_set(idx_snr, :) = true_velocity_ref_rad_row;

    % 仅用于脚本末尾显示最后一次有效估计及其对应真实值，不参与 RMSE 统计。
    last_valid_trial = find(~isnan(nor_delay_esti_row), 1, 'last');
    if ~isempty(last_valid_trial)
        nor_delay_esti = nor_delay_esti_row(last_valid_trial);
        nor_Doppler_esti = nor_Doppler_esti_row(last_valid_trial);
        P_sig_strong_esti = P_sig_strong_esti_row(last_valid_trial);
        range_esti_rad = range_esti_rad_row(last_valid_trial);
        velocity_esti_rad = velocity_esti_rad_row(last_valid_trial);
        true_nor_delay = true_nor_delay_ref_row(last_valid_trial);
        true_nor_Doppler = true_nor_Doppler_ref_row(last_valid_trial);
        true_range_ref_rad = true_range_ref_rad_row(last_valid_trial);
        true_velocity_ref_rad = true_velocity_ref_rad_row(last_valid_trial);
        idx_ref_target = ref_target_index_row(last_valid_trial);
    end

    if cfg.enableDiagnosticPlots && idx_snr == cfg.diagnosticSNRIndex
        diagnostic_ML_value_refine = diagnostic_ML_value_row{cfg.diagnosticTrialIndex};
        diagnostic_parameters_esti = diagnostic_parameters_row{cfg.diagnosticTrialIndex};
        if isempty(diagnostic_ML_value_refine) || isempty(diagnostic_parameters_esti)
            warning('Diagnostic trial SNR index %d trial index %d did not produce a valid estimate.', ...
                cfg.diagnosticSNRIndex, cfg.diagnosticTrialIndex);
        else
            diagnostic_nor_delay_esti = diagnostic_parameters_esti(1);
            diagnostic_nor_Doppler_esti = diagnostic_parameters_esti(2);
            diagnostic_P_sig_strong_esti = diagnostic_parameters_esti(3);

            %% ===================== 反演强目标回波并计算残差 =====================
            L_delay_noDelayDop = 2;
            L_Doppler_noDelayDop = 2;
            L_delay_aver = 1 / (Nor_delay_Search_set(2) - Nor_delay_Search_set(1));
            L_Doppler_aver = 1 / (Nor_Doppler_Search_set(2) - Nor_Doppler_Search_set(1));

            nor_delay_esti_strong_upSam = diagnostic_nor_delay_esti * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
            nor_Doppler_esti_strong_upSam = diagnostic_nor_Doppler_esti * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);

            D_func_search_tmp = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor( ...
                AmbigFunc_pulsedAFDM_noDelayDop, para.N, para.L, L_delay_noDelayDop, L_Doppler_noDelayDop, ...
                diagnostic_P_sig_strong_esti, nor_delay_esti_strong_upSam, nor_Doppler_esti_strong_upSam);

            % 从高分辨率反演结果中抽取与当前接收搜索网格一致的区域
            echo_strong_esti = D_func_search_tmp(1:L_delay_noDelayDop:max_norDelay*para.L*L_delay_noDelayDop, ...
                para.N*para.L*L_Doppler_noDelayDop/2 - cfg.max_norDoppler*para.L*L_Doppler_noDelayDop/2 + 1 : ...
                para.N*para.L*L_Doppler_noDelayDop/2 + cfg.max_norDoppler*para.L*L_Doppler_noDelayDop/2);

            ML_value_refine_set_aver_weak = abs(diagnostic_ML_value_refine - echo_strong_esti);

            [value_row_max_refine, I_row_set_max_refine] = max(diagnostic_ML_value_refine);
            [~, I_col_max_refine] = max(value_row_max_refine);
            I_row_max_refine = I_row_set_max_refine(I_col_max_refine);

            % 补偿徙动画图
            figure;
            [X, Y] = meshgrid(Nor_Doppler_Search_set, Nor_delay_Search_set);
            mesh(X, Y, 10 * log10(diagnostic_ML_value_refine / max(max(diagnostic_ML_value_refine))));
            ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
            xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
            zlabel('Normalized average squared MF');
            set(gca, 'FontSize', 10);
            axis([-inf inf -inf inf -40 2]);
            title('单目标-补偿徙动-实际匹配滤波结果-非相干累积');

            figure;
            mesh(X, Y, 10 * log10(abs(ML_value_refine_set_aver_weak) ./ max(max(abs(echo_strong_esti)))));
            ylabel('Normalized delay \tau');
            xlabel('Normalized Doppler \nu');
            axis([-inf inf -inf inf -40 0]);
            title('Residual after strong target cancellation');

            figure;
            plot(Nor_delay_Search_set, 10 * log10(abs(diagnostic_ML_value_refine(:, I_col_max_refine)) ./ ...
                max(max(abs(echo_strong_esti)))), '-b*');
            hold on;
            plot(Nor_delay_Search_set, 10 * log10(abs(echo_strong_esti(:, I_col_max_refine)) ./ ...
                max(max(abs(echo_strong_esti)))), '--ro');
            grid on;
            axis([-inf inf -30 2]);
            set(gca, 'FontSize', 10);
            xlabel('Normalized delay \tau');
            ylabel('Normalized average squared MF');
            legend('Received', 'Reconstructed', 'Location', 'south');
            annotation('textbox', [0.45 0.8 0.1 0.1], 'String', ['hat tau = ', num2str(diagnostic_nor_delay_esti, '%10.5f')], ...
                'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10);
            title('Delay slice comparison');

            figure;
            plot(Nor_Doppler_Search_set, 10 * log10(abs(diagnostic_ML_value_refine(I_row_max_refine, :)) ./ ...
                max(max(abs(echo_strong_esti)))), '-b*');
            hold on;
            plot(Nor_Doppler_Search_set, 10 * log10(abs(echo_strong_esti(I_row_max_refine, :)) ./ ...
                max(max(abs(echo_strong_esti)))), '--ro');
            grid on;
            axis([-inf inf -30 2]);
            set(gca, 'FontSize', 10);
            xlabel('Normalized Doppler \nu');
            ylabel('Normalized average squared MF');
            legend('Received', 'Reconstructed', 'Location', 'south');
            annotation('textbox', [0.45 0.8 0.1 0.1], 'String', ['hat nu = ', num2str(diagnostic_nor_Doppler_esti, '%10.5f')], ...
                'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10);
            title('Doppler slice comparison');
        end
    end
    fprintf('Finished SNR %+g dB: valid trials = %d/%d, failed trials = %d/%d\n', ...
        SNR_dB_curr, sum(~isnan(range_error_m_row)), num_trial_run, sum(detect_failed_row), num_trial_run);
end

%% ===================== 4. 蒙特卡洛 统计、画图与保存 =====================
nor_delay_rmse = nan(num_snr, 1);
nor_Doppler_rmse = nan(num_snr, 1);
range_rmse_m = nan(num_snr, 1);
velocity_rmse_mps = nan(num_snr, 1);

% RMSE 只统计成功产生估计的 trial；检测失败率由 detect_failed 单独统计。
for idx_snr = 1:num_snr
    valid_nor_delay = ~isnan(nor_delay_error(idx_snr, :));
    valid_nor_Doppler = ~isnan(nor_Doppler_error(idx_snr, :));
    valid_range = ~isnan(range_error_m(idx_snr, :));
    valid_velocity = ~isnan(velocity_error_mps(idx_snr, :));
    if any(valid_nor_delay)
        nor_delay_rmse(idx_snr) = sqrt(mean(nor_delay_error(idx_snr, valid_nor_delay).^2));
    end
    if any(valid_nor_Doppler)
        nor_Doppler_rmse(idx_snr) = sqrt(mean(nor_Doppler_error(idx_snr, valid_nor_Doppler).^2));
    end
    if any(valid_range)
        range_rmse_m(idx_snr) = sqrt(mean(range_error_m(idx_snr, valid_range).^2));
    end
    if any(valid_velocity)
        velocity_rmse_mps(idx_snr) = sqrt(mean(velocity_error_mps(idx_snr, valid_velocity).^2));
    end
end
detect_fail_rate = mean(detect_failed, 2);

% mc_result 同时保存 trial 级原始估计、误差矩阵和 SNR 级 RMSE/失败率统计。
mc_result = struct();
mc_result.SNR_dB_set = SNR_run_set;
mc_result.numMonteCarlo = num_trial_run;
mc_result.true_nor_delay_all = true_nor_delay_set;
mc_result.true_nor_Doppler_all = true_nor_Doppler_set;
mc_result.true_range_rad = true_range_rad_set;
mc_result.true_velocity_rad = true_velocity_rad_set;
mc_result.ref_target_index = ref_target_index_set;
mc_result.true_nor_delay_ref = true_nor_delay_ref_set;
mc_result.true_nor_Doppler_ref = true_nor_Doppler_ref_set;
mc_result.true_range_ref_rad = true_range_ref_rad_set;
mc_result.true_velocity_ref_rad = true_velocity_ref_rad_set;
mc_result.nor_delay_esti = nor_delay_esti_set;
mc_result.nor_Doppler_esti = nor_Doppler_esti_set;
mc_result.P_sig_strong_esti = P_sig_strong_esti_set;
mc_result.range_esti_rad = range_esti_rad_set;
mc_result.velocity_esti_rad = velocity_esti_rad_set;
mc_result.nor_delay_error = nor_delay_error;
mc_result.nor_Doppler_error = nor_Doppler_error;
mc_result.range_error_m = range_error_m;
mc_result.velocity_error_mps = velocity_error_mps;
mc_result.nor_delay_rmse = nor_delay_rmse;
mc_result.nor_Doppler_rmse = nor_Doppler_rmse;
mc_result.range_rmse_m = range_rmse_m;
mc_result.velocity_rmse_mps = velocity_rmse_mps;
mc_result.detect_failed = detect_failed;
mc_result.detect_fail_rate = detect_fail_rate;
mc_result.trial_message = trial_message;

fprintf('\n========== 蒙特卡洛 RMSE Summary ==========\n');
if isnan(idx_ref_target)
    fprintf('No valid trial estimate was produced; reference target summary is unavailable.\n');
else
    fprintf('Last valid reference target index = %d / %d\n', idx_ref_target, cfg.Num_target);
    fprintf('Last valid reference true delay   = %.6f\n', true_nor_delay);
    fprintf('Last valid reference true Doppler = %.6f\n', true_nor_Doppler);
    fprintf('Last valid reference true range   = %.6f m\n', true_range_ref_rad);
    fprintf('Last valid reference true velocity = %.6f m/s\n', true_velocity_ref_rad);
end
for idx_snr = 1:num_snr
    fprintf(['SNR = %+g dB: norm-delay RMSE = %.6g, norm-Doppler RMSE = %.6g, ' ...
        'range RMSE = %.6g m, velocity RMSE = %.6g m/s, fail rate = %.2f%%\n'], ...
        SNR_run_set(idx_snr), nor_delay_rmse(idx_snr), nor_Doppler_rmse(idx_snr), ...
        range_rmse_m(idx_snr), velocity_rmse_mps(idx_snr), 100 * detect_fail_rate(idx_snr));
end
fprintf('==============================================\n\n');

if cfg.enableSummaryPlots
    figure;
    semilogy(SNR_run_set, nor_delay_rmse, '-o', 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)');
    ylabel('Normalized Delay RMSE');
    title('Normalized Delay Estimation RMSE vs SNR');

    figure;
    semilogy(SNR_run_set, nor_Doppler_rmse, '-o', 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)');
    ylabel('Normalized Doppler RMSE');
    title('Normalized Doppler Estimation RMSE vs SNR');

    figure;
    semilogy(SNR_run_set, range_rmse_m, '-o', 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)');
    ylabel('Range RMSE (m)');
    title('Range Estimation RMSE vs SNR');

    figure;
    semilogy(SNR_run_set, velocity_rmse_mps, '-o', 'LineWidth', 1.5);
    grid on;
    xlabel('SNR (dB)');
    ylabel('Velocity RMSE (m/s)');
    title('Velocity Estimation RMSE vs SNR');
end

if cfg.saveMonteCarloResult
    result_dir = fileparts(cfg.monteCarloResultPath);
    if ~isempty(result_dir) && ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    save(cfg.monteCarloResultPath, 'cfg', 'para', 'mc_result');
    fprintf('蒙特卡洛 result saved to %s\n', cfg.monteCarloResultPath);
end

fprintf('AFDM 蒙特卡洛 SNR sweep finished.\n');
fprintf('Last estimated normalized delay   = %.6f\n', nor_delay_esti);
fprintf('Last estimated normalized Doppler = %.6f\n', nor_Doppler_esti);
toc;
