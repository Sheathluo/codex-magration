% function [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v12(factor_delay_refine, factor_Doppler_refine, ...
%     Tx_signal, Rx_signal, Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler, L, Nc)
function [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v12(Tx_signal, Rx_signal, Nor_delay_set, Nor_Doppler_set, L, Nc)
%%  版本信息
%% 2026/09/15 V12.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026  BUAA
%======================================================
% Designer：Ni Yuanhan
% Function Description：在时域利用最大似然估计算法，快速估计目标时延、多普勒
% Input parameter:  Tx_signal - 发送波形
%                   Rx_signal - 接收波形
%                   Min_delay - 需要搜索的最小时延
%                   Max_delay - 需要搜索的时延范围[0,Max_delay]
%                   Min_nor_Doppler - 需要搜索的最小多普勒
%                   Max_nor_Doppler - 需要搜索的时延范围[-Max_Doppler,Max_Doppler]
%                   factor_delay_refine - 时延搜索间隔
%                   factor_Doppler_refine - 多普勒搜索间隔
%                   L - 过采样倍数
%                   Nc - 子载波数
% Output parameter: ML_value_refine - 二维匹配滤波功率面，[grid_delay * grid_dopple]
% Modification: 移除了未使用的变量，利用矩阵乘法彻底取代了双重for循环，极大提升执行效率。
%               在v10的基础上把搜索网格的生成移到函数外面。 2026-09-15，姓名：LuoTao
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1. 生成搜索网格
    % Nor_delay_set = Min_delay : factor_delay_refine : Max_delay;
    % Nor_Doppler_set = Min_nor_Doppler : factor_Doppler_refine : Max_nor_Doppler;

    Len_delay = length(Nor_delay_set);
    % Len_Doppler = length(Nor_Doppler_set);
    NL = Nc * L;

    % 2. 预计算时延维度矩阵 (Rx_intDelay_mat)
    % 大小: [NL x Len_delay]
    Rx_intDelay_mat = zeros(NL, Len_delay);
    for idx_delay = 1:Len_delay
        nor_delay = Nor_delay_set(idx_delay);
        Rx_intDelay_mat(:, idx_delay) = Rx_signal(1 + nor_delay : NL + nor_delay, 1);
    end

    % 3. 预计算多普勒维度矩阵 (Rx_copy_mat)
    % 大小: [NL x Len_Doppler]
    % 利用矩阵运算一次性生成所有 Doppler 相位补偿
    time_vec = (0 : NL - 1).';
    Doppler_phase = exp(-1i * 2 * pi / NL * (time_vec * Nor_Doppler_set));

    % 利用隐式扩展，将 Tx_signal 的共轭应用到整个多普勒相位矩阵上
    Rx_copy_mat = conj(Tx_signal) .* Doppler_phase;

    % 4. 矩阵乘法计算最终的最大似然统计量面 (ML_value_refine)
    % 巧妙利用转置矩阵乘积，完美等效并替代原先的双重 for 循环中的点乘求和
    % 大小: [Len_delay x Len_Doppler]
    ML_value_refine = (1 / NL) * (Rx_intDelay_mat.' * Rx_copy_mat);

end