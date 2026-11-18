function [delay_grid, doppler_grid] = f_RFT_searchRange_v1(SNR_dB, true_delay, ...
    true_nor_doppler, Nsym, Min_delay_full, Max_delay_full, Min_nor_doppler_full, ...
    Max_nor_doppler_full, factor_delay_refine, factor_doppler_refine)
%% 版本信息
%% 2026/09/14 V1.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：根据 SNR、累计符号数和真实 delay-Doppler 参数，自适应
%              设置 f_AFDM_RFT_v3 的搜索中心、搜索范围和搜索步长。
%
% Modification:
%   f_AFDM_RFT_searchRange_v1 V1.0: - 由真实参数确定栅格中心，由 SNR 和累计符号数确定搜索范围；
%   f_AFDM_RFT_searchRange_v1 V1.1: - 删除诊断信息输出，并增加时延和Doppler 搜索步长的自适应调整；
%   f_AFDM_RFT_searchRange_v1 V1.2: - 将时延和 Doppler 最大搜索步长分别限制为 1 和 0.5；
%   f_AFDM_RFT_searchRange_v1 V1.3: - 直接输出时延和 Doppler 搜索网格，供 f_AFDM_RFT_v3 内部调用；
%   f_AFDM_RFT_searchRange_v1 V1.4: - 低于 SNR 阈值时采用完整搜索范围；局部搜索固定点数，只调整步长；
%   f_AFDM_RFT_searchRange_v1 V1.5: - 修正局部搜索的时延步长，确保最大为 1；
%                                     修改日期：2026年09月14日 修改人：Codex
%
% 输入参数：
%   SNR_dB                       - 接收信号 SNR，单位 dB
%   true_delay                   - 真实归一化时延，单位为上采样采样点
%   true_nor_Doppler             - 真实归一化 Doppler，nu = fd/delta_f
%   Nsym                         - RFT 相干累计的 AFDM 符号数
%   Min_delay_full, Max_delay_full - 完整时延搜索范围
%   Min_nor_Doppler_full, Max_nor_doppler_full - 完整 Doppler 搜索范围
%   factor_delay_refine, factor_doppler_refine - 完整的时延和多普勒搜索步长
%
% 输出参数：
%   delay_grid                   - RFT 时延搜索网格
%   doppler_grid                 - RFT 归一化 Doppler 搜索网格
%
% 说明：
%   本函数使用真实参数构造 oracle 栅格，仅用于仿真加速。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 以 500 个累计符号为参考；累计符号数每增加 4 倍，指标提高约 6 dB。
    SNR_metric_dB = SNR_dB + 10 * log10(Nsym / 500);

    % 输入 SNR 低于该阈值时，不使用真实值缩小范围，直接搜索完整网格。
    SNR_full_search_threshold_dB = -30;
    if SNR_dB <= SNR_full_search_threshold_dB
        delay_grid = Min_delay_full:factor_delay_refine:Max_delay_full;
        doppler_grid = Min_nor_doppler_full:factor_doppler_refine:Max_nor_doppler_full;
        return;
    end

    % 局部搜索始终保持 11 个时延点和 21 个 Doppler 点。
    delay_half_grid_count = 50;
    Doppler_half_grid_count = 50;

    % 每行依次为：[SNR上界, delay步长, Doppler步长, 保留小数位]；时延步长最大为 1。
    search_table = [...
        -25,    0.01,        0.01,   2, 2; ...
        -20,    0.01,       0.01,   2, 2; ...
        -15,    0.01,       0.01,   2, 2; ...
        -10,    0.01,      0.01,   2, 2; ...
        -5 ,    0.01,      0.01,   2, 2; ...
        inf,    0.01,    0.01,     2, 2;];

    table_index = find(SNR_metric_dB <= search_table(:, 1), 1, 'first');
    factor_delay = search_table(table_index, 2);
    factor_Doppler = search_table(table_index, 3);
    decimal_delay = search_table(table_index, 4);
    decimal_doppler = search_table(table_index, 5);

    % 以真实参数为中心生成固定点数网格，越界时整体平移而不减少点数。
    delay_center_round = round(true_delay, decimal_delay);
    delay_start = delay_center_round - delay_half_grid_count * factor_delay;
    delay_end = delay_start + 2 * delay_half_grid_count * factor_delay;
    if delay_start < Min_delay_full
        delay_start = delay_start + ceil((Min_delay_full - delay_start) / factor_delay) * factor_delay;
    elseif delay_end > Max_delay_full
        delay_start = delay_start - ceil((delay_end - Max_delay_full) / factor_delay) * factor_delay;
    end
    delay_grid = delay_start + (0:2 * delay_half_grid_count) * factor_delay;

    true_nor_doppler_round = round(true_nor_doppler, decimal_doppler);
    Doppler_start = true_nor_doppler_round - Doppler_half_grid_count * factor_Doppler;
    Doppler_end = Doppler_start + 2 * Doppler_half_grid_count * factor_Doppler;
    if Doppler_start < Min_nor_doppler_full
        Doppler_start = Doppler_start + ...
            ceil((Min_nor_doppler_full - Doppler_start) / factor_Doppler - 1e-12) * factor_Doppler;
    elseif Doppler_end > Max_nor_doppler_full
        Doppler_start = Doppler_start - ...
            ceil((Doppler_end - Max_nor_doppler_full) / factor_Doppler - 1e-12) * factor_Doppler;
    end
    doppler_grid = Doppler_start + (0:2 * Doppler_half_grid_count) * factor_Doppler;
    % doppler_grid_round = round(doppler_grid, decimal_doppler);
end
