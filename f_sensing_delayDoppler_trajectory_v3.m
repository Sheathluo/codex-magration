function [tau_vector] = f_sensing_delayDoppler_trajectory_v3(cfg, para, samples_per_symbol, doppler_search)
%%  版本信息
%%
%%2026/06/16 V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：生成距离徙动补偿用的时延轨迹，并计算真实距离和速度
%
% 输入参数：
%   cfg                - Monte Carlo 和目标配置结构体
%   cfg.Num_target     - 真实目标/路径数量
%   cfg.chan_coef      - 各目标复增益，长度应等于 cfg.Num_target
%   cfg.doppler_i_cal  - 各目标归一化多普勒，长度应等于 cfg.Num_target
%   cfg.delay_i_cal    - 各目标初始归一化时延，长度应等于 cfg.Num_target
%   cfg.flag_delayMigration - 1 表示符号级徙动量，2 表示采样点级徙动量，3 表示零徙动
%   para               - AFDM 和雷达物理参数结构体
%   samples_per_symbol - 单个 pulse-shaped AFDM 符号的基带采样点数
%
% 输出参数：
%   tau_vector      - 各目标逐采样点归一化“徙动量”轨迹，尺寸为 cfg.Num_target x signalLen
%                     注意：该输出不包含初始时延 delay_i_cal，适合做徙动补偿
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 由归一化多普勒换算成目标速度，用于计算随时间增长的归一化时延漂移
    true_velocity = doppler_search * para.delta_f / para.fc * para.c / 2;

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
        tau = repmat(doppler_search, 1, para.Nsym) + delta_tau_migration;
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

    % 这里用于去保护区后的接收矩阵补偿，因此不再额外拼接 preSuffix/head/tail。
    tau_vector_head = repmat(tau(:, 1), 1, para.N_preSuffix);
    tau_vector_end = repmat(tau(:, end), 1, para.N_preSuffix);
    tau_vector = [tau_vector_head, tau_vector_sym, tau_vector_end];
    tau_vector = tau;
end
