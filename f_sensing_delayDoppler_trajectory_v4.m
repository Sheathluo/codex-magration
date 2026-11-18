function [chan_coef, doppler_i_cal, delay_i_cal, tau_vector, delay_max, true_range, true_velocity] = ...
    f_sensing_delayDoppler_trajectory_v4(cfg, para, samples_per_symbol)
%%  版本信息
%%
%%2026/07/14 V4.0
% Modify:在增加符号级徙动时，避免使用符号开头或结尾作为整个符号都时延，使用符号中间点的时延作为整个符号都时延
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：生成感知目标时延/多普勒轨迹，并计算真实距离和速度
%
% 输入参数：
%   cfg                - Monte Carlo 和目标配置结构体
%   cfg.Num_target     - 真实目标/路径数量
%   cfg.chan_coef      - 各目标复增益，长度应等于 cfg.Num_target
%   cfg.doppler_i_cal  - 各目标归一化多普勒，长度应等于 cfg.Num_target
%   cfg.delay_i_cal    - 各目标初始归一化时延，长度应等于 cfg.Num_target
%   cfg.flag_delayMigration - 1 表示符号级距离徙动，2 表示采样点级距离徙动
%   para               - AFDM 和雷达物理参数结构体
%   samples_per_symbol - 单个 pulse-shaped AFDM 符号的基带采样点数
%
% 输出参数：
%   chan_coef       - cfg.chan_coef 的列向量形式，各目标复增益
%   v_i_cal         - cfg.v_i_cal 的列向量形式，各目标归一化多普勒
%   delay_i_cal     - cfg.delay_i_cal 的列向量形式，各目标初始归一化时延
%   tau_vector      - 各目标逐采样点归一化时延轨迹，尺寸为 cfg.Num_target x signalLen
%   delay_max       - 信道输出需要预留的最大整数时延采样数
%   true_range      - 各目标真实双程雷达距离，单位 m，尺寸为 cfg.Num_target x 1
%   true_velocity   - 各目标真实双程雷达速度，单位 m/s，尺寸为 cfg.Num_target x 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 统一转成列向量，后续按目标维度做矩阵化距离徙动轨迹。
    chan_coef = cfg.chan_coef(:);
    doppler_i_cal = cfg.doppler_vec(:);
    delay_i_cal = cfg.delay_vec(:);

    % if length(chan_coef) ~= cfg.Num_target || length(doppler_i_cal) ~= cfg.Num_target || length(delay_i_cal) ~= cfg.Num_target
    %     error('cfg.Num_target must match chan_coef, v_i_cal, and delay_i_cal lengths.');
    % end

    % 根据CP确定最大时延
    max_initial_delay = (para.Ncp - 1) * para.L - 1;
    % if any(delay_i_cal < 0) || any(delay_i_cal > max_initial_delay)
    %     error('All initial target delays must be within the receiver search range [0, (para.Ncp - 1) * para.L - 1].');
    % end

    % 由归一化时延换算成距离
    true_range = delay_i_cal / para.fs_up * para.c / 2;
    % 由归一化多普勒换算成目标速度，用于计算随时间增长的归一化时延漂移
    true_velocity = doppler_i_cal * para.delta_f / para.fc * para.c / 2;

    if cfg.flag_delayMigration == 1
        % 符号级距离徙动：每个 AFDM 符号内时延保持不变，符号间更新。
        % delta_tau_migration = true_velocity * 2 * ((1:para.Nsym) * para.L * samples_per_symbol / para.c);  % 第一个符号加了徙动，不合理，所以注释
        % delta_tau_migration = true_velocity * 2 * ((0:para.Nsym-1) * para.L * samples_per_symbol / para.c);  % 用了符号的第一个采样点的时延代替整个符号的时延，不合理
        % 使用符号中间的采样点作为符号的时延
        %%%% 时间向量
        % time_vec = ((0:para.Nsym-1) * para.L * samples_per_symbol + floor(para.L * samples_per_symbol / 2));
        time_vec = (0:para.Nsym-1) * para.L * samples_per_symbol; %%%% !! 暂时使用符号的第一个点作为整个符号的时延
        %%%% 双程距离变化
        delta_range = true_velocity * 2 * time_vec;
        %%%% 时延变化
        delta_tau_migration = delta_range / para.c;
        tau = repmat(delay_i_cal, 1, para.Nsym) + delta_tau_migration;
        tau_vector_sym = kron(tau, ones(1, para.L * samples_per_symbol));
    elseif cfg.flag_delayMigration == 2
        % 采样点级距离徙动：每个采样点都有对应的时延轨迹。
        % delta_tau_migration = true_velocity * ((1:(para.Nsym * para.L * samples_per_symbol)) * 2 / para.c);
        %%%% 时间向量
        time_vec = 0:(para.Nsym * para.L * samples_per_symbol)-1;
        %%%% 双程距离变化
        delta_range = true_velocity * 2 * time_vec;
        %%%% 时延变化
        delta_tau_migration = delta_range / para.c;
        tau = repmat(delay_i_cal, 1, para.Nsym * para.L * samples_per_symbol) + delta_tau_migration;
        tau_vector_sym = tau;
    elseif cfg.flag_delayMigration == 3
        % 不产生徙动
        true_velocity = 0;
        % delta_tau_migration = true_velocity * ((1:para.Nsym) * 2 * para.L * samples_per_symbol / para.c);
        delta_tau_migration = true_velocity * ((0:para.Nsym-1) * 2 * para.L * samples_per_symbol / para.c);
        tau = repmat(delay_i_cal, 1, para.Nsym) + delta_tau_migration;
        tau_vector_sym = kron(tau, ones(1, para.L * samples_per_symbol));
    else
        error('cfg.flag_delayMigration 必须为 1（符号级）或 2（采样点级）或 3（不徙动）。');
    end

    % 前后保护区沿用首尾时延，保证 tau_vector 与带保护区的发射信号长度对齐。
    tau_vector_head = repmat(tau(:, 1), 1, para.N_preSuffix);
    tau_vector_end = repmat(tau(:, end), 1, para.N_preSuffix);
    tau_vector = [tau_vector_head, tau_vector_sym, tau_vector_end];
    % delay_max = para.Ncp * para.L - 1;
    delay_max = para.Ncp * para.L;

end
