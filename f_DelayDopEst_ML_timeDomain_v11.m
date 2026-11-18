function [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v11(Tx_signal, Rx_signal, nor_delay_grid, nor_doppler, para)
%%  版本信息
%% 2026/09/14 V11.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026  BUAA
%======================================================
% Designer：Luo Tao
% Function Description：二维匹配滤波
% Input parameter:  Tx_signal - 发送波形
%                   Rx_signal - 接收目标回波
%                   delay_grid - 需要搜索的时延网格，向量
%                   nor_doppler - 需要搜索的多普勒值，标量
%                   para - 波形参数结构体
% Output parameter: ML_value_refine - 二维匹配滤波功率面，[delay_grid * 1(nor_doppler)]
% Modification:
%  ---------------------------------------------------------------------------------
%  Version      Time        Modifier               Modification Content
%  ---------------------------------------------------------------------------------
%  v10       --             [罗韬]     [移除了未使用的变量，利用矩阵乘法彻底取代了双重for循环，极大提升执行效率。]
%  v11    2026-09-14        [罗韬]     [修改内容:v11版在v10版基础上：
%                                       1、精简了输出输出参数；
%                                       2、利用分数时延滤波器实现参考信号的分数时延，以支持分数时延搜索]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1. 生成搜索网格
    % Nor_delay_set = Min_delay : factor_delay_refine : Max_delay;
    % Nor_Doppler_set = Min_nor_Doppler : factor_Doppler_refine : Max_nor_Doppler;

    len_delay = length(nor_delay_grid);
    NL = para.L * para.N;

    % 2. 预计算时延维度矩阵 (Rx_intDelay_mat), 大小: [NL x Len_delay]
    Rx_fracDelay_mat = zeros(NL, len_delay);
    % 提取整数时延
    nor_int_delays = floor(nor_delay_grid);
    % 获取分数时延
    nor_fracd_delays = nor_delay_grid - nor_int_delays;
    % 把不同符号的时延扩展为矩阵
    nor_fracd_delays_mat = kron(nor_fracd_delays, ones(NL, 1));

    for delay_idx = 1:len_delay
        nor_int_delay = nor_int_delays(delay_idx);
        % 实现整数时延
        Rx_intDelay = Rx_signal(1 + nor_int_delay : NL + nor_int_delay, 1);
        nor_frac_delay = nor_fracd_delays_mat(:, delay_idx);
        % 实现分数时延
        Rx_fracDelay_mat(:,delay_idx) = f_tv_delay_intfrac_v3(Rx_intDelay, -nor_frac_delay, para.sincFilterLen);
    end

    % 3. 预计算多普勒维度矩阵 (Rx_copy_mat)
    % 大小: [NL x Len_Doppler]
    % 利用矩阵运算一次性生成所有 Doppler 相位补偿
    time_vec = (0 : NL - 1).';
    Doppler_phase = exp(-1i * 2 * pi / NL * (time_vec * nor_doppler));

    % 利用隐式扩展，将 Tx_signal 的共轭应用到整个多普勒相位矩阵上
    Rx_copy_mat = conj(Tx_signal) .* Doppler_phase;

    % 4. 矩阵乘法计算最终的最大似然统计量面 (ML_value_refine)
    % 巧妙利用转置矩阵乘积，完美等效并替代原先的双重 for 循环中的点乘求和
    % 大小: [Len_delay x Len_Doppler]
    ML_value_refine = (1 / NL) * (Rx_fracDelay_mat.' * Rx_copy_mat);

end