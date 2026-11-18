function [g_ps, P_rrc] = f_pulse_shaping_filter(para)
%%  版本信息
%%
%%2026/06/16 V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：生成 AFDM 脉冲成形滤波器和接收端等效脉冲成形矩阵
%
% 输入参数：
%   para          - AFDM 与脉冲成形参数结构体
%   para.flag_RRC - 1 表示使用根升余弦滤波器；否则使用升余弦滤波器
%   para.beta    - 滚降因子
%   para.L       - 过采样倍数
%   para.Nc      - 子载波数/有效 AFDM 符号长度
%   para.N_preSuffix - 每个 AFDM 符号前后保护采样长度
%   para.N_tapFiler  - 脉冲成形滤波器抽头跨度参数
%
% 输出参数：
%   g_ps  - 归一化脉冲成形滤波器，作用于上采样后的带 CP/保护区发射序列
%   P_rrc - 等效循环脉冲成形矩阵，用于生成无 CP 有效 AFDM 符号的匹配滤波参考
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 生成发射端使用的单位范数 RRC 脉冲成形滤波器。
    g_ps = rcosdesign(para.beta, para.N_tapFiler, para.L, 'sqrt').';
    g_ps = g_ps / norm(g_ps);
    % % % 将脉冲嵌入一个有效 AFDM 符号长度内，形成循环移位矩阵的第一列。
    p_rrc_trunc = zeros(para.Nc * para.L, 1);
    start_idx_rrc = para.Nc * para.L / 2 - para.N_preSuffix * para.L / 2 + 1;
    end_idx_rrc = start_idx_rrc + length(g_ps) - 1;

    p_rrc_trunc(start_idx_rrc:end_idx_rrc) = g_ps;
    p_rrc_ = circshift(p_rrc_trunc, -para.Nc * para.L / 2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % p_rrc_ambig = rcosdesign(para.beta, para.Nc, para.L, 'sqrt').';
    % p_rrc_ambig = p_rrc_ambig(1:end-1);
    %
    % p_rrc_ambig_trunc = zeros(para.Nc * para.L, 1);
    % idx_start = para.Nc * para.L / 2 - para.N_preSuffix * para.L + 1;
    % idx_end = idx_start + (2 * para.N_preSuffix * para.L + 1) - 1;
    % p_rrc_ambig_trunc(idx_start:idx_end) = p_rrc_ambig(idx_start:idx_end);
    % p_rrc_ambig = circshift(p_rrc_ambig_trunc, -para.Nc * para.L / 2);
    % p_rrc_ = p_rrc_ambig;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % P_rrc 的每一列是 p_rrc 的循环移位，对应一个有效采样点的脉冲响应。
    P_rrc = zeros(para.Nc * para.L, para.Nc * para.L);
    for idx = 1:para.Nc * para.L
        P_rrc(:, idx) = circshift(p_rrc_, idx - 1);
    end
end
