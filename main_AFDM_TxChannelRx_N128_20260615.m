%% AFDM 发射 - 时变感知信道 - 接收处理一体化主链路
% 纯脚本版本：不在主文件末尾写本地函数，方便后续直接改参数、插断点和分段调试。
% 原三段脚本中的发射/参考/接收 dat 文件交互已改为内存变量直接传递。

clear;
close all;
rng(1);
tic;

%% ===================== 0. 参数配置 =====================
cfg = struct();
cfg.enablePlots = true;        % 是否画图。批量跑性能时建议改为 false
cfg.enableNoise = false;       % 是否给信道输出加噪。默认保持原信道脚本中“生成噪声但不加噪”的行为
cfg.useParallel = true;        % 是否使用 parfor 做逐符号匹配滤波
cfg.rngSeed = 1;               % 固定随机种子，保证发射数据可复现

% 感知信道目标参数。这里默认沿用 main_sensingChannel_N128_20260511 的单目标设置
cfg.Num_target = 1;            % 仿真的真实目标数/路径数
cfg.Max_num_target_est = 2;    % 接收端 OMP 最多检测目标数，不一定等于真实目标数
cfg.chan_coef = 1;             % 目标/路径复增益
cfg.v_i_cal = 11.72;           % 归一化多普勒
cfg.delay_i_cal = 10.35;       % 归一化时延，支持非整数时延
cfg.SNR_in_dB = 10;            % 若 cfg.enableNoise=true，则按该 SNR 加噪
cfg.sincFilterLen = 500;       % 时变分数时延 sinc 插值滤波器长度
cfg.flag_delayMigration = 2;   % 1: 符号级距离徙动；2: 采样点级距离徙动

% 接收端截取与搜索参数
cfg.index_start = 0;           % 接收信号起始截取点
cfg.index_start_noise = 1e3;   % 噪声估计起始偏移，仿真中通常只作为备用估计
cfg.rx_L_delay = 1;            % 接收端时延搜索步进为 1/rx_L_delay
cfg.max_norDoppler = 15;       % 接收端粗搜索最大归一化多普勒

% 理论模糊函数计算较慢，因此仍允许缓存。该缓存不是发射-信道-接收的中间交互文件
cfg.ambigCachePath = './result/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v6_20251022.mat';

para = struct();
% 脉冲成形参数
para.L = 2;                    % 过采样倍数
para.beta = 0.35;              % RRC/RC 滚降因子
para.L_delay = 2;              % 理论模糊函数的时延维过采样倍数
para.L_Doppler = 2;            % 理论模糊函数/接收搜索的多普勒维过采样倍数
para.flag_pulse = 1;           % 1 表示使用脉冲成形
para.flag_RRC = 1;             % 1 表示使用根升余弦滤波器
para.N_preSuffix = 20;         % 每个 AFDM 符号前后额外保护点长度
para.N_tapFiler = 20;          % 成形滤波器半边抽头参数
para.fs_up = 80e6 / 8;         % 上采样后采样率
para.fs = para.fs_up / para.L;

% AFDM/雷达基本参数
para.c = 3e8;
para.fc = 24e9;
para.Nc = 128;                 % 子载波数
para.N = para.Nc;
para.M = para.Nc;
para.Nsym = 500;               % 符号数/帧数
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

L = para.L;
Nc = para.Nc;
N = para.N;
Nsym = para.Nsym;
Ncp = para.Ncp;
N_preSuffix = para.N_preSuffix;
samples_per_symbol = Nc + Ncp + 2 * N_preSuffix;

%% ===================== 1. 发射端：生成 AFDM 波形 =====================
rng(cfg.rngSeed);

% 生成随机比特并做 16-QAM 调制。这里沿用原发射脚本的数据组织方式
bit_num = Nc;
data_bits = randi([0 1], bit_num * Nsym * log2(para.Mod), 1);
mod_para = struct('Mod', para.Mod, 'Tb', 1/20/1e6, 'fs', 20e6, 'Rb', 20e6);

data_bits_reshape = reshape(data_bits, log2(para.Mod), []).';
data_symbols = bi2de(data_bits_reshape);
qam_symbols = qammod(data_symbols, para.Mod, 'UnitAveragePower', true).';
kron_vec = ones(1, round(mod_para.fs * mod_para.Tb));
qam_upsamp = kron(qam_symbols, kron_vec);
qam_frame = reshape(qam_upsamp, length(qam_upsamp) / Nsym, []);

waveform_index = 'AFDM';
[~, transmit_data_noCP_AFDM, ~, ~] = f_waveform_generator_ADDM_v3( ...
    waveform_index, para.c1, para.c2, para.M, para.N, Ncp, Nsym, qam_frame);

% 按原发射脚本“先加额外前后缀和 CP，再过采样，再脉冲成形”的方式构造发射波形
transmit_data_CP_preSuffix_AFDM = [ ...
    transmit_data_noCP_AFDM(end-Ncp-N_preSuffix+1:end-Ncp, :); ...
    transmit_data_noCP_AFDM(end-Ncp+1:end, :); ...
    transmit_data_noCP_AFDM; ...
    transmit_data_noCP_AFDM(1:N_preSuffix, :)];

x_sigCP_preSuffix = reshape(transmit_data_CP_preSuffix_AFDM, [], 1).';
x_upSamCP_preSuffix = reshape([sqrt(L) * x_sigCP_preSuffix; ...
    zeros(L-1, samples_per_symbol * Nsym)], [], 1);

if para.flag_pulse == 1
    if para.flag_RRC == 1
        g_ps = rcosdesign(para.beta, para.N_tapFiler, L, 'sqrt').';
        g_ps = g_ps / norm(g_ps);
    else
        g_ps = rcosdesign(para.beta, 2 * para.N_tapFiler, L, 'normal').';
        g_ps = g_ps / norm(g_ps);
    end
else
    g_ps = [zeros(para.N_tapFiler * L, 1); 1; zeros(para.N_tapFiler * L, 1)];
end

x_pulseShaped_CP_preSuffix = conv(g_ps, x_upSamCP_preSuffix);
tx_signal = 0.1 * x_pulseShaped_CP_preSuffix ./ max(abs(x_pulseShaped_CP_preSuffix));

% 构造接收端匹配滤波使用的等效脉冲成形参考信号。原脚本中该信号曾写成 Tx_ConfencedSignal dat 文件
p_rrc_trunc = zeros(Nc * L, 1);
start_idx_rrc = Nc * L / 2 - N_preSuffix * L / 2 + 1;
end_idx_rrc = start_idx_rrc + length(g_ps) - 1;
if end_idx_rrc > length(p_rrc_trunc)
    error('Reference pulse exceeds one AFDM symbol. Check N_preSuffix or N_tapFiler.');
end
p_rrc_trunc(start_idx_rrc:end_idx_rrc) = g_ps;
p_rrc = circshift(p_rrc_trunc, -Nc * L / 2);

P_rrc = zeros(Nc * L, Nc * L);
for idx = 1:Nc * L
    P_rrc(:, idx) = circshift(p_rrc, idx - 1);
end

x_upSam = kron(sqrt(L) * transmit_data_noCP_AFDM, [1; zeros(L-1, 1)]);
Tx_signal_pulsedAFDM_eff = P_rrc * x_upSam;

if cfg.enablePlots
    figure;
    plot(20 * log10(abs(tx_signal)));
    title('Tx signal spectrum-like magnitude');

    figure;
    plot(abs(tx_signal));
    title('Tx signal envelope');
end

%% ===================== 2. 感知信道：非整数时延 + 多普勒 + 距离徙动 =====================
chan_coef = cfg.chan_coef(:);
v_i_cal = cfg.v_i_cal(:);
delay_i_cal = cfg.delay_i_cal(:);

if length(chan_coef) ~= cfg.Num_target || length(v_i_cal) ~= cfg.Num_target || length(delay_i_cal) ~= cfg.Num_target
    error('cfg.Num_target must match chan_coef, v_i_cal, and delay_i_cal lengths.');
end

velocity = v_i_cal .* para.delta_f * para.c / para.fc;

% 距离徙动建模方式：
% 模式1对应原脚本中的 delta_tau = (1:Nsym) * 2 .* (velocity' .* L * (Nc + Ncp + 2 * N_preSuffix) / c);
%       即符号间时延变化，符号内时延保持不变。
% 模式2对应原脚本中注释掉的逐采样点公式：
%       delta_tau = (1:(Nsym * L * (Nc + Ncp + 2 * N_preSuffix))) * 2 .* (velocity' / c);
%       即符号内每个采样点的时延也随速度连续变化。
if cfg.flag_delayMigration == 1
    delta_tau_migration = velocity * ((1:Nsym) * 2 * L * samples_per_symbol / para.c);
    tau = repmat(delay_i_cal, 1, Nsym) + delta_tau_migration;
    tau_vector_sym = kron(tau, ones(1, L * samples_per_symbol));
elseif cfg.flag_delayMigration == 2
    delta_tau_migration = velocity * ((1:(Nsym * L * samples_per_symbol)) * 2 / para.c);
    tau = repmat(delay_i_cal, 1, Nsym * L * samples_per_symbol) + delta_tau_migration;
    tau_vector_sym = tau;
else
    error('cfg.flag_delayMigration 必须为 1（符号级）或 2（采样点级）。');
end

% f_time_vary_channel_v3 按采样点接收时延轨迹，因此统一补齐前后保护区的 tau_vector
tau_vector_head = repmat(tau(:, 1), 1, N_preSuffix);
tau_vector_end = repmat(tau(:, end), 1, N_preSuffix);
tau_vector = [tau_vector_head, tau_vector_sym, tau_vector_end];
delay_max = Ncp * L - 1;

[Rx_signal_AFDM_tmp, ~] = f_time_vary_channel_v3(tx_signal(:), Nc, L, cfg.Num_target, ...
    chan_coef, tau_vector, v_i_cal, delay_max, cfg.sincFilterLen, 1);

P_rx_cal = mean(abs(Rx_signal_AFDM_tmp).^2);
P_noise = P_rx_cal / (10^(cfg.SNR_in_dB / 10));
len_rx_channel = length(Rx_signal_AFDM_tmp);
if para.flag_RRC == 1
    noise_AFDM = sqrt(L) * sqrt(P_noise / 2) * (randn(len_rx_channel, 1) + 1i * randn(len_rx_channel, 1));
else
    noise_AFDM = sqrt(P_noise / 2) * (randn(len_rx_channel, 1) + 1i * randn(len_rx_channel, 1));
end

if cfg.enableNoise
    Rx_signal_AFDM_tmp = Rx_signal_AFDM_tmp + noise_AFDM;
end

%% ===================== 3. 接收端初始化与理论模糊函数 =====================
rx_L_delay = cfg.rx_L_delay;
factor_delay_refine = 1 / rx_L_delay;
factor_Doppler_refine = 1 / para.L_Doppler;

% 搜索区域：时延只搜 CP 范围内，多普勒按 cfg.max_norDoppler 限定
max_norDelay = Ncp - 1;
max_norDoppler = cfg.max_norDoppler;
Min_delay = 0;
Max_delay = max_norDelay * L - 1;
Min_nor_Doppler = -max_norDoppler * L / 2;
Max_nor_Doppler = max_norDoppler * L / 2 - 1 / para.L_Doppler;

Nor_delay_refineSearch_set = Min_delay:factor_delay_refine:Max_delay;
Nor_Doppler_refineSearch_set = Min_nor_Doppler:factor_Doppler_refine:Max_nor_Doppler;
len_delay = length(Nor_delay_refineSearch_set);
len_doppler = length(Nor_Doppler_refineSearch_set);

cache_dir = fileparts(cfg.ambigCachePath);
if ~isempty(cache_dir) && ~exist(cache_dir, 'dir')
    mkdir(cache_dir);
end

if exist(cfg.ambigCachePath, 'file')
    ambig_data = load(cfg.ambigCachePath, 'AmbigFunc_pulsedAFDM_noDelayDop');
    AmbigFunc_pulsedAFDM_noDelayDop = ambig_data.AmbigFunc_pulsedAFDM_noDelayDop;
else
    % 缓存不存在时现场生成理论模糊函数。该步骤可能较慢
    if para.flag_pulse == 1
        p_rrc_ambig = rcosdesign(para.beta, Nc, L, 'sqrt').';
        p_rrc_ambig = p_rrc_ambig(1:end-1);
    else
        p_rrc_ambig = [1; zeros(Nc * L - 1, 1)];
    end

    p_rrc_ambig_trunc = zeros(Nc * L, 1);
    idx_start = Nc * L / 2 - N_preSuffix * L + 1;
    idx_end = idx_start + (2 * N_preSuffix * L + 1) - 1;
    p_rrc_ambig_trunc(idx_start:idx_end) = p_rrc_ambig(idx_start:idx_end);
    p_rrc_ambig = circshift(p_rrc_ambig_trunc, -Nc * L / 2);

    delay_nor_set = 0:1/rx_L_delay:N * L - 1/rx_L_delay;
    Doppler_nor_set = -N * L / 2:1/para.L_Doppler:N * L / 2 - 1/para.L_Doppler;

    AmbigFunc_pulsedAFDM_noDelayDop = f_theoretical_AmbigFunc_pulseShaped_AFDM_v7( ...
        p_rrc_ambig, N, L, para.c1, para.mu_4, delay_nor_set, Doppler_nor_set, 0, 0);
    save(cfg.ambigCachePath, 'AmbigFunc_pulsedAFDM_noDelayDop');
end
AmbigFunc_pulsedAFDM_noDelayDop = AmbigFunc_pulsedAFDM_noDelayDop / (N^2);

%% ===================== 4. 接收端预处理：截取、去过渡带、移除保护间隔 =====================
samples_per_symbol_L = samples_per_symbol * L;
Len_sig = samples_per_symbol_L * Nsym + N_preSuffix * L;

if cfg.index_start + Len_sig > length(Rx_signal_AFDM_tmp)
    error('Requested Rx segment exceeds received signal length. Check Nsym or cfg.index_start.');
end
Rx_signal_AFDM = Rx_signal_AFDM_tmp(cfg.index_start + 1 : cfg.index_start + Len_sig);

% 去掉成形滤波器前后过渡带
start_trunc = para.N_tapFiler * L / 2 + 1;
end_trunc = para.N_tapFiler * L / 2 + samples_per_symbol_L * Nsym;
Rx_signal_AFDM_trunc = Rx_signal_AFDM(start_trunc:end_trunc);

% reshape 后一次性按符号取出有效区域，替代原脚本逐符号 for 循环截取
rx_frame = reshape(Rx_signal_AFDM_trunc, samples_per_symbol_L, Nsym);
valid_start = (Ncp + N_preSuffix) * L + 1;
valid_idx = valid_start : valid_start + (N + Ncp) * L - 1;
X_pulseShaped_remove_preSuffix = rx_frame(valid_idx, :);

if cfg.index_start_noise + Len_sig <= length(Rx_signal_AFDM_tmp)
    noise_received = Rx_signal_AFDM_tmp(cfg.index_start_noise + Len_sig : end);
    P_noise_cal = mean(abs(noise_received).^2);
else
    n_noise = min(100, length(Rx_signal_AFDM_tmp));
    P_noise_cal = mean(abs(Rx_signal_AFDM_tmp(1:n_noise)).^2);
end

%% ===================== 5. 粗距离徙动补偿 + 时域匹配滤波 =====================
Nor_Doppler_set = Nor_Doppler_refineSearch_set;
est_doppler_len = length(Nor_Doppler_set);
ML_value_refine_set_vactor = zeros(len_delay, len_doppler, est_doppler_len);

disp('Running coarse range-migration compensation and time-domain matched filtering...');
for dopple_index = 1:est_doppler_len
    % 对每个候选多普勒，先按估计速度粗补偿距离徙动，再做逐符号匹配滤波
    velocity_search = Nor_Doppler_set(dopple_index) .* para.delta_f * para.c / para.fc;
    migration_delta_tau_norm = (1:Nsym) * 2 .* (velocity_search .* L * samples_per_symbol / para.c);
    migration_delta_tau_time = migration_delta_tau_norm / para.fs;

    if Nor_Doppler_set(dopple_index) >= 0
        range_migration_bin = floor(migration_delta_tau_time / para.delta_tau);
    else
        range_migration_bin = ceil(migration_delta_tau_time / para.delta_tau);
    end

    X_circ = zeros(size(X_pulseShaped_remove_preSuffix));
    for n_sym = 1:Nsym
        X_circ(:, n_sym) = circshift(X_pulseShaped_remove_preSuffix(:, n_sym), -range_migration_bin(n_sym));
    end

    ML_value_refine_set_curr = zeros(len_delay, len_doppler, Nsym);
    if cfg.useParallel
        parfor n_sym = 1:Nsym
            ML_value_refine = f_DelayDopEst_ML_timeDomain_v10(factor_delay_refine, factor_Doppler_refine, ...
                Tx_signal_pulsedAFDM_eff(:, n_sym), X_circ(:, n_sym), ...
                Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler, L, N);
            ML_value_refine_set_curr(:, :, n_sym) = ML_value_refine;
        end
    else
        for n_sym = 1:Nsym
            ML_value_refine = f_DelayDopEst_ML_timeDomain_v10(factor_delay_refine, factor_Doppler_refine, ...
                Tx_signal_pulsedAFDM_eff(:, n_sym), X_circ(:, n_sym), ...
                Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler, L, N);
            ML_value_refine_set_curr(:, :, n_sym) = ML_value_refine;
        end
    end

    ML_value_refine_set_vactor(:, :, dopple_index) = sum(abs(ML_value_refine_set_curr).^2, 3) ./ Nsym;
end

% 选择非相干累积峰值最大的粗多普勒页，作为后续精估计输入
[~, max_idx] = max(ML_value_refine_set_vactor(:));
[~, ~, max_page] = ind2sub(size(ML_value_refine_set_vactor), max_idx);
ML_value_refine_set_aver = ML_value_refine_set_vactor(:, :, max_page);
[X, Y] = meshgrid(Nor_Doppler_refineSearch_set, Nor_delay_refineSearch_set);

if cfg.enablePlots
    figure;
    mesh(X, Y, 10 * log10(ML_value_refine_set_aver / max(max(ML_value_refine_set_aver))));
    ylabel('Normalized delay \tau');
    xlabel('Normalized Doppler \nu');
    zlabel('Normalized average squared MF');
    set(gca, 'FontSize', 10);
    axis([-inf inf -inf inf -40 2]);
    title('Matched filter after coarse range-migration compensation');
end

%% ===================== 6. 梯度下降/OMP 精细估计非整数时延和多普勒 =====================
% 理论模糊函数反折，等价于原脚本中的 X(-tau,-nu)，这里用索引直接实现
AmbigFunc_flip = AmbigFunc_pulsedAFDM_noDelayDop([1, end:-1:2], [1, end:-1:2]);

N_delay_peak = 4;
N_Doppler_peak = 4;
Num_iter_single = 1;
Num_iter_cyclic = 30;
Pfa_target = 1e-3;
P_noise_MF = P_noise_cal ./ N;
Th_noTarget = 1 / Nsym * gaminv(1 - Pfa_target, Nsym, P_noise_MF);

disp('Running fine parameter estimation with theoretical ambiguity function...');
set_nor_parameters_esti = f_paraEstibyGradDesc_newtonOMP_multiTarget_v1( ...
    ML_value_refine_set_aver, AmbigFunc_flip, Nor_delay_refineSearch_set, Nor_Doppler_refineSearch_set, ...
    N, L, rx_L_delay, para.L_Doppler, max_norDelay, max_norDoppler, N_delay_peak, N_Doppler_peak, ...
    para.fs_up, para.delta_f, Num_iter_single, Num_iter_cyclic, Th_noTarget, cfg.Max_num_target_est);

if isempty(set_nor_parameters_esti)
    error('No target was detected. Check SNR, target parameters, and detection threshold.');
end

nor_delay_esti = set_nor_parameters_esti(1, 1);
nor_Doppler_esti = set_nor_parameters_esti(2, 1);
P_sig_strong_esti = set_nor_parameters_esti(3, 1);
range_esti_rad = nor_delay_esti / para.fs_up * para.c / 2;
velocity_esti_rad = nor_Doppler_esti * para.delta_f / para.fc * para.c / 2;

%% ===================== 7. 反演强目标回波并计算残差 =====================
L_delay_noDelayDop = 2;
L_Doppler_noDelayDop = 2;
L_delay_aver = 1 / (Nor_delay_refineSearch_set(2) - Nor_delay_refineSearch_set(1));
L_Doppler_aver = 1 / (Nor_Doppler_refineSearch_set(2) - Nor_Doppler_refineSearch_set(1));

nor_delay_esti_strong_upSam = nor_delay_esti * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
nor_Doppler_esti_strong_upSam = nor_Doppler_esti * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);

D_func_search_tmp = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor( ...
    AmbigFunc_pulsedAFDM_noDelayDop, N, L, L_delay_noDelayDop, L_Doppler_noDelayDop, ...
    P_sig_strong_esti, nor_delay_esti_strong_upSam, nor_Doppler_esti_strong_upSam);

% 从高分辨率反演结果中抽取与当前接收搜索网格一致的区域
echo_strong_esti = D_func_search_tmp(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop, ...
    N*L*L_Doppler_noDelayDop/2 - max_norDoppler*L*L_Doppler_noDelayDop/2 + 1 : ...
    N*L*L_Doppler_noDelayDop/2 + max_norDoppler*L*L_Doppler_noDelayDop/2);

ML_value_refine_set_aver_weak = abs(ML_value_refine_set_aver - echo_strong_esti);

[value_row_max_refine, I_row_set_max_refine] = max(ML_value_refine_set_aver);
[~, I_col_max_refine] = max(value_row_max_refine);
I_row_max_refine = I_row_set_max_refine(I_col_max_refine);

%% ===================== 8. 绘图与对外输出变量 =====================
if cfg.enablePlots
    figure;
    mesh(X, Y, 10 * log10(abs(echo_strong_esti) ./ max(max(abs(echo_strong_esti)))));
    ylabel('Normalized delay \tau');
    xlabel('Normalized Doppler \nu');
    title('Reconstructed strong target echo');

    figure;
    mesh(X, Y, 10 * log10(abs(ML_value_refine_set_aver_weak) ./ max(max(abs(echo_strong_esti)))));
    ylabel('Normalized delay \tau');
    xlabel('Normalized Doppler \nu');
    axis([-inf inf -inf inf -40 0]);
    title('Residual after strong target cancellation');

    figure;
    plot(Nor_delay_refineSearch_set, 10 * log10(abs(ML_value_refine_set_aver(:, I_col_max_refine)) ./ ...
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
    annotation('textbox', [0.45 0.8 0.1 0.1], 'String', ['hat tau = ', num2str(nor_delay_esti, '%10.5f')], ...
        'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10);
    title('Delay slice comparison');

    figure;
    plot(Nor_Doppler_refineSearch_set, 10 * log10(abs(ML_value_refine_set_aver(I_row_max_refine, :)) ./ ...
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
    annotation('textbox', [0.45 0.8 0.1 0.1], 'String', ['hat nu = ', num2str(nor_Doppler_esti, '%10.5f')], ...
        'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10);
    title('Doppler slice comparison');
end

fprintf('AFDM integrated chain finished.\n');
fprintf('Estimated normalized delay   = %.6f\n', nor_delay_esti);
fprintf('Estimated normalized Doppler = %.6f\n', nor_Doppler_esti);
toc;
