%% AFDM 蒙特卡洛 SNR 扫描：发射 - 时变感知信道 - 接收估计
% 主脚本保留实验配置、Monte Carlo 主循环、统计、画图和保存。
% 发射、接收预处理、匹配滤波等链路模块由 f_* 函数承接。

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
cfg.enablePlots = false;                        % 保留兼容字段；parfor 内不画图，诊断图由 enableDiagnosticPlots 控制
cfg.addNoiseFlag = false;                        % 是否加噪
% cfg.SNR_dB_set = -10:5:30;                      % 蒙特卡洛扫描的 SNR 点
cfg.SNR_dB_set = 10;                           % 蒙特卡洛扫描的 SNR 点
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

rng('shuffle');                % 蒙特卡洛仿真不固定随机种子
% rng(1);
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

% 感知信道目标参数
cfg.Max_num_target_est = 2;    % 接收端 OMP 最多检测目标数，不一定等于真实目标数
cfg.chan_coef = [1];             % 目标/路径复增益

% 根据归一化值计算速度
decimals = 2;                   % 时延和多普勒保留小数位数
cfg.v_i_cal = [0.3];           % 归一化多普勒
cfg.delay_i_cal = [10.64];       % 归一化时延，支持非整数时延
cfg.range = cfg.delay_i_cal./para.fs*para.c/2;                % 距离，m
cfg.velocity = cfg.v_i_cal.*para.delta_f*para.c/(2*para.fc);  % 速度，m/s

% 根据速度计算归一化值
% cfg.range = 310.5;             % 距离，m
% cfg.velocity = 2861.32;        % 速度，m/s
% cfg.delay_i_cal = cfg.range .* 2 / para.c * para.fs;
% cfg.v_i_cal = cfg.velocity .* 2 * para.fc / para.c / para.delta_f;

cfg.sincFilterLen = 300;       % 时变分数时延 sinc 插值滤波器长度
cfg.flag_delayMigration = 3;   % 1: 符号级距离徙动；2: 采样点级距离徙动；3：不徙动
fprintf('\n================ 徙动模式 ===================\n');
fprintf('距离徙动模式 flag_delayMigration       = %s', delayMigrationStr(cfg.flag_delayMigration));
fprintf('\n============================================\n');

cfg.chan_coef = cfg.chan_coef(:);
cfg.v_i_cal = cfg.v_i_cal(:);
cfg.delay_i_cal = cfg.delay_i_cal(:);
if length(cfg.chan_coef) ~= length(cfg.delay_i_cal) || length(cfg.v_i_cal) ~= length(cfg.delay_i_cal)
    error('cfg.chan_coef, cfg.v_i_cal, and cfg.delay_i_cal must have the same length.');
end
cfg.Num_target = length(cfg.delay_i_cal);

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

% 生成目标的时变时延轨迹 tau_vector，并返回真实距离/速度作为蒙特卡洛误差基准。
[chan_coef, v_i_cal, delay_i_cal, tau_vector, delay_max, true_range_rad, true_velocity_rad] = ...
    f_sensing_delayDoppler_trajectory_v1(cfg, para, samples_per_symbol);

% 多目标场景下只统计真实最强目标；估计端仍使用第一列检测结果作为对比对象。
[~, idx_ref_target] = max(abs(chan_coef));
true_nor_delay = delay_i_cal(idx_ref_target);
true_nor_Doppler = v_i_cal(idx_ref_target);
true_range_ref_rad = true_range_rad(idx_ref_target);
true_velocity_ref_rad = true_velocity_rad(idx_ref_target);

factor_delay_refine = 1 / cfg.rx_L_delay;
factor_Doppler_refine = 1 / para.L_Doppler;

max_norDelay = para.Ncp - 1;
Min_delay = 0;
Max_delay = max_norDelay * para.L - 1;
Min_nor_Doppler = -cfg.max_norDoppler * para.L / 2;
Max_nor_Doppler = cfg.max_norDoppler * para.L / 2 - 1 / para.L_Doppler;

Nor_delay_refineSearch_set = Min_delay:factor_delay_refine:Max_delay;
Nor_Doppler_refineSearch_set = Min_nor_Doppler:factor_Doppler_refine:Max_nor_Doppler;

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
range_error_m = nan(num_snr, num_trial_run);
velocity_error_mps = nan(num_snr, num_trial_run);
detect_failed = false(num_snr, num_trial_run);
trial_message = cell(num_snr, num_trial_run);

nor_delay_esti = NaN;
nor_Doppler_esti = NaN;
P_sig_strong_esti = NaN;
range_esti_rad = NaN;
velocity_esti_rad = NaN;

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
    range_error_m_row = nan(1, num_trial_run);
    velocity_error_mps_row = nan(1, num_trial_run);
    detect_failed_row = false(1, num_trial_run);
    trial_message_row = cell(1, num_trial_run);
    diagnostic_ML_value_row = cell(1, num_trial_run);
    diagnostic_parameters_row = cell(1, num_trial_run);

    % parfor idx_trial = 1:num_trial_run
    for idx_trial = 1:num_trial_run
        trial_timer = tic;
        try
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
            % ML_value_refine_set_aver = f_AFDM_coarseRM_MF_v1(X_pulseShaped_remove_preSuffix, ...
            %     Tx_signal_pulsedAFDM_eff, para, Nor_Doppler_refineSearch_set, factor_delay_refine, ...
            %     factor_Doppler_refine, Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler);

            ML_value_refine_set_aver = f_AFDM_coarseRM_MF_v2(X_pulseShaped_remove_preSuffix, ...
                Tx_signal_pulsedAFDM_eff, para, Nor_Doppler_refineSearch_set, factor_delay_refine, ...
                factor_Doppler_refine, Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler);

            % 用 v4 返回的噪声功率设置检测门限，使门限和实际加噪功率口径一致。
            P_noise_MF = P_noise_cal ./ para.N;
            Th_noTarget = 1 / para.Nsym * gaminv(1 - Pfa_target, para.Nsym, P_noise_MF);

            % 基于匹配滤波面和理论模糊函数做 OMP/梯度精估计，输出归一化时延、多普勒和强目标功率。
            set_nor_parameters_esti = f_paraEstibyGradDesc_newtonOMP_multiTarget_v2( ...
                ML_value_refine_set_aver, AmbigFunc_flip, Nor_delay_refineSearch_set, Nor_Doppler_refineSearch_set, ...
                para.N, para.L, cfg.rx_L_delay, para.L_Doppler, max_norDelay, cfg.max_norDoppler, N_delay_peak, N_Doppler_peak, ...
                para.fs_up, para.delta_f, Num_iter_single, Num_iter_cyclic, Th_noTarget, cfg.Max_num_target_est);

            if isempty(set_nor_parameters_esti)
                detect_failed_row(idx_trial) = true;
                trial_message_row{idx_trial} = 'No target was detected.';
                send(progress_queue, sprintf('[no target] SNR %+g dB | trial %d/%d | elapsed %.1f s', ...
                    SNR_dB_curr, idx_trial, num_trial_run, toc(trial_timer)));
            else
                nor_delay_esti_curr = round(set_nor_parameters_esti(1, 1),decimals);
                nor_Doppler_esti_curr = round(set_nor_parameters_esti(2, 1),decimals);
                P_sig_strong_esti_curr = round(set_nor_parameters_esti(3, 1),decimals);
                range_esti_rad_curr = nor_delay_esti_curr / para.fs_up * para.c / 2;
                velocity_esti_rad_curr = nor_Doppler_esti_curr * para.delta_f / para.fc * para.c / 2;

                % 记录当前 trial 的原始估计值和相对真实目标参数的误差
                nor_delay_esti_row(idx_trial) = nor_delay_esti_curr;
                nor_Doppler_esti_row(idx_trial) = nor_Doppler_esti_curr;
                P_sig_strong_esti_row(idx_trial) = P_sig_strong_esti_curr;
                range_esti_rad_row(idx_trial) = range_esti_rad_curr;
                velocity_esti_rad_row(idx_trial) = velocity_esti_rad_curr;
                range_error_m_row(idx_trial) = range_esti_rad_curr - true_range_ref_rad;
                velocity_error_mps_row(idx_trial) = velocity_esti_rad_curr - true_velocity_ref_rad;
                send(progress_queue, sprintf(['[done] SNR %+g dB | trial %d/%d | delay %.4f | Doppler %.4f ' ...
                    '| range %.3f m | velocity %.3f m/s | elapsed %.1f s'], ...
                    SNR_dB_curr, idx_trial, num_trial_run, nor_delay_esti_curr, nor_Doppler_esti_curr, ...
                    range_esti_rad_curr, velocity_esti_rad_curr, toc(trial_timer)));

                if cfg.enableDiagnosticPlots && idx_snr == cfg.diagnosticSNRIndex && idx_trial == cfg.diagnosticTrialIndex
                    diagnostic_ML_value_row{idx_trial} = ML_value_refine_set_aver;
                    diagnostic_parameters_row{idx_trial} = [nor_delay_esti_curr; nor_Doppler_esti_curr; P_sig_strong_esti_curr];
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
    range_error_m(idx_snr, :) = range_error_m_row;
    velocity_error_mps(idx_snr, :) = velocity_error_mps_row;
    detect_failed(idx_snr, :) = detect_failed_row;
    trial_message(idx_snr, :) = trial_message_row;

    last_valid_trial = find(~isnan(nor_delay_esti_row), 1, 'last');
    if ~isempty(last_valid_trial)
        nor_delay_esti = nor_delay_esti_row(last_valid_trial);
        nor_Doppler_esti = nor_Doppler_esti_row(last_valid_trial);
        P_sig_strong_esti = P_sig_strong_esti_row(last_valid_trial);
        range_esti_rad = range_esti_rad_row(last_valid_trial);
        velocity_esti_rad = velocity_esti_rad_row(last_valid_trial);
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
            L_delay_aver = 1 / (Nor_delay_refineSearch_set(2) - Nor_delay_refineSearch_set(1));
            L_Doppler_aver = 1 / (Nor_Doppler_refineSearch_set(2) - Nor_Doppler_refineSearch_set(1));

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
            [X, Y] = meshgrid(Nor_Doppler_refineSearch_set, Nor_delay_refineSearch_set);
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
            plot(Nor_delay_refineSearch_set, 10 * log10(abs(diagnostic_ML_value_refine(:, I_col_max_refine)) ./ ...
                max(max(abs(echo_strong_esti)))), '-b*');
            hold on;
            plot(Nor_delay_refineSearch_set, 10 * log10(abs(echo_strong_esti(:, I_col_max_refine)) ./ ...
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
            plot(Nor_Doppler_refineSearch_set, 10 * log10(abs(diagnostic_ML_value_refine(I_row_max_refine, :)) ./ ...
                max(max(abs(echo_strong_esti)))), '-b*');
            hold on;
            plot(Nor_Doppler_refineSearch_set, 10 * log10(abs(echo_strong_esti(I_row_max_refine, :)) ./ ...
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
range_rmse_m = nan(num_snr, 1);
velocity_rmse_mps = nan(num_snr, 1);

% RMSE 只统计成功产生估计的 trial；检测失败率由 detect_failed 单独统计。
for idx_snr = 1:num_snr
    valid_range = ~isnan(range_error_m(idx_snr, :));
    valid_velocity = ~isnan(velocity_error_mps(idx_snr, :));
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
mc_result.true_nor_delay_all = delay_i_cal;
mc_result.true_nor_Doppler_all = v_i_cal;
mc_result.true_nor_delay = true_nor_delay;
mc_result.true_nor_Doppler = true_nor_Doppler;
mc_result.true_range_rad = true_range_rad;
mc_result.true_velocity_rad = true_velocity_rad;
mc_result.ref_target_index = idx_ref_target;
mc_result.true_nor_delay_ref = true_nor_delay;
mc_result.true_nor_Doppler_ref = true_nor_Doppler;
mc_result.true_range_ref_rad = true_range_ref_rad;
mc_result.true_velocity_ref_rad = true_velocity_ref_rad;
mc_result.nor_delay_esti = nor_delay_esti_set;
mc_result.nor_Doppler_esti = nor_Doppler_esti_set;
mc_result.P_sig_strong_esti = P_sig_strong_esti_set;
mc_result.range_esti_rad = range_esti_rad_set;
mc_result.velocity_esti_rad = velocity_esti_rad_set;
mc_result.range_error_m = range_error_m;
mc_result.velocity_error_mps = velocity_error_mps;
mc_result.range_rmse_m = range_rmse_m;
mc_result.velocity_rmse_mps = velocity_rmse_mps;
mc_result.detect_failed = detect_failed;
mc_result.detect_fail_rate = detect_fail_rate;
mc_result.trial_message = trial_message;

fprintf('\n========== 蒙特卡洛 RMSE Summary ==========\n');
fprintf('Reference target index = %d / %d\n', idx_ref_target, cfg.Num_target);
fprintf('Reference true range    = %.6f m\n', true_range_ref_rad);
fprintf('Reference true velocity = %.6f m/s\n', true_velocity_ref_rad);
for idx_snr = 1:num_snr
    fprintf('SNR = %+g dB: range RMSE = %.6g m, velocity RMSE = %.6g m/s, fail rate = %.2f%%\n', ...
        SNR_run_set(idx_snr), range_rmse_m(idx_snr), velocity_rmse_mps(idx_snr), 100 * detect_fail_rate(idx_snr));
end
fprintf('==============================================\n\n');

if cfg.enableSummaryPlots
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
    save(cfg.monteCarloResultPath, 'cfg', 'para', 'mc_result', 'true_range_rad', 'true_velocity_rad');
    fprintf('蒙特卡洛 result saved to %s\n', cfg.monteCarloResultPath);
end

fprintf('AFDM 蒙特卡洛 SNR sweep finished.\n');
fprintf('Last estimated normalized delay   = %.6f\n', nor_delay_esti);
fprintf('Last estimated normalized Doppler = %.6f\n', nor_Doppler_esti);
toc;
