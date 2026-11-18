function E_mat_total = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v5( ...
    X, P, norDelay_strong, norDopper_strong, para, echo_rows, echo_cols, cfg_trial)
%% 版本信息
% Copyright (C) 2025 BUAA
% 原始函数设计：Ni Yuanhan；v2/v3 修改：LuoTao。
% 2026/10/23 V4.0：按相同补偿位移分组，先累加频域相位再生成理论模板。
% 2026/11/06：公共多普勒移位后先抽取所需列，减少各组时延 FFT 的规模。
% 2026/11/10 V5.0：在组内相位求和时加入逐符号主峰功率权重。
% 权重由实际发射参考能量平方计算，主脚本归一化为平均值1。
% 保留v4分组、先截取再移位及除以符号数的处理顺序。
%
% 输入：
%   X                 - 离线生成的无时延、无多普勒理论模糊函数矩阵
%   P                 - 目标功率
%   norDelay_strong   - 每个符号在理论模板网格上的时延，行/列向量均可
%   norDopper_strong  - 所有符号共用的多普勒，单位为理论模板网格点
%   para.Nsym        - 积累符号数 K
%   echo_rows/cols    - 从完整理论模板抽取接收搜索区域的索引
%   cfg_trial         - flag_delayMigration、2-by-K migration_compensate_bin
%                       以及逐符号权重 symbol_power_weight（平均值为1）
% 输出：
%   E_mat_total       - 抽取、补偿并跨符号平均后的实值二维响应
%
% 原理：本函数生成加权期望功率面，不重新计算统计模糊函数。
% 每个符号使用固定权重w_k：A = sum_k w_k * S_k(A_k) / K。
% 权重是功率权重，直接乘w_k，不取平方根。
% 时延向和多普勒向移位作用于不同维度，可以交换顺序。
% 多普勒移位对所有符号相同，因此只计算一次。对具有相同整数补偿
% 位移的符号，抽取和 circshift 也是相同的线性算子，可先累加时延
% 相位，再统一逆 FFT、抽取和移位。不同补偿位移之间仍分别计算，
% 不把接收搜索区域上的 circshift 合并到完整模板的频域相位中。
%
% 设模板有 S=M_tau*M_nu 个点、G 种补偿位移、C 个接收多普勒列：
% 主要计算量由 v3 的 O(K*S*log(S)) 降为
% O(S*log(M_nu)+K*M_tau+(G+1)*M_tau*C*log(M_tau))。
% G 接近 K 时分组收益变小；不对任意输入承诺固定加速倍数。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% 1. 准备逐符号时延和补偿位移分组
    num_symbol = para.Nsym;
    symbol_weights = cfg_trial.symbol_power_weight(1:num_symbol);
    symbol_weights = symbol_weights(:);
    no_migration = cfg_trial.flag_delayMigration == 3;
    delay_values = norDelay_strong(:);
    % 无徙动模式仅使用第一个时延；仍保留各符号的补偿算子。
    if no_migration
        delay_values = repmat(delay_values(1), num_symbol, 1);
    end
    % v3 只读取前 K 个符号；额外的时延和补偿列不参与当前帧计算。
    delay_values = delay_values(1:num_symbol);
    compensation_shift = -cfg_trial.migration_compensate_bin(1:2, 1:num_symbol).';
    [unique_shift, ~, group_index] = unique(compensation_shift, 'rows', 'stable');
    num_group = size(unique_shift, 1);
    group_members = accumarray(group_index, (1:num_symbol).', [num_group, 1], @(indices) {indices});

    %% 2. 公共多普勒移位只计算一次，再缓存时延向频谱
    [num_delay, num_doppler] = size(X);
    delay_bins_centered = ((0:num_delay - 1) - num_delay / 2).';
    doppler_bins_centered = (0:num_doppler - 1) - num_doppler / 2;
    doppler_phase = exp(-1i * 2 * pi / num_doppler * norDopper_strong * doppler_bins_centered);

    % 保留 v3 的两次 fftshift 约定，避免改变其频率 bin 和移位定义。
    doppler_spectrum_centered = fftshift(fft(P * X, [], 2), 2);
    template_doppler_shifted = ifft(fftshift(doppler_spectrum_centered .* doppler_phase, 2), [], 2);
    % 多普勒移位仍在完整模板上完成。后续只沿时延维运算，因此可提前
    % 抽取所需多普勒列；列顺序及重复索引保持不变，尚不做 circshift。
    template_doppler_selected = template_doppler_shifted(:, echo_cols);
    delay_spectrum_centered = fftshift(fft(template_doppler_selected, [], 1), 1);

    %% 3. 组内先累加相位，每组仅生成一次理论模板
    response_acc = zeros(numel(echo_rows), numel(echo_cols));
    if no_migration
        common_delay_phase = exp(-1i * 2 * pi / num_delay * delay_values(1) * delay_bins_centered);
    end
    for group_idx = 1:num_group
        members = group_members{group_idx};
        if no_migration
            delay_phase_sum = sum(symbol_weights(members)) * common_delay_phase;
        else
            % 非共轭转置：每列对应一个符号，沿符号维累加。
            delay_phase_group = exp((-1i * 2 * pi / num_delay * delay_bins_centered) * delay_values(members).');
            delay_phase_sum = sum(delay_phase_group .* symbol_weights(members).', 2);
        end
        response_group = ifft(fftshift(delay_spectrum_centered .* delay_phase_sum, 1), [], 1);

        % 必须先抽取，再在接收搜索网格上移位，保持与 v3 完全相同的算子。
        response_search = response_group(echo_rows, :);
        response_search = circshift(response_search, unique_shift(group_idx, :));
        response_acc = response_acc + response_search;
    end

    %% 4. 跨符号平均，舍弃浮点运算留下的虚部
    E_mat_total = real(response_acc / num_symbol);
end
