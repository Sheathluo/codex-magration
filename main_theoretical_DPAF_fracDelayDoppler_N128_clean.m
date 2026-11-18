%% 生成理论 DPAF（非整数时延 / 非整数多普勒）—— 整理版
% 来源：main_generate_theoretical_DPAF_fracDelayDoppler_N128_20251124.m
% 设计：Ni Yuanhan (BUAA, 2025)   整理：2026-11-12
% 功能：生成 RRC 脉冲成形 AFDM 在无时延、无多普勒回波下的理论模糊函数。
% 注意：源脚本计算结果仅保留在工作区，save 语句处于注释状态（见文末，取消注释即可保存）。
%
% 依赖（同目录）：
%   f_theoretical_AmbigFunc_pulseShaped_AFDM_v6.m  —— 自定义函数，已提供
%   rcosdesign  —— MATLAB 内置（Signal Processing Toolbox）
%   parfor      —— 函数内部使用（Parallel Computing Toolbox）
%   （目录中另有 v7 向量化提速版，接口与 v6 完全相同，如需切换只改一行函数名）
%
% 相对源文件删除的内容（均为未被使用、不影响结果的代码或变量）：
%   1) 全部注释掉的旧版本代码（v2/v3 调用、旧 .mat 的 save/load、绘图代码）
%   2) 死变量：L_delay, L_Doppler, flag_RRC, N_tapFiler（定义未用，截断实际用
%      N_preSuffix）, Nsym, Mod, M, N_delay, N_Doppler, tic（无对应 toc）
%   3) c1/c2 的被覆盖中间赋值，生效值仅 c1 = 8/(2*Nc)；c2 未传入函数，整块删除
%      （含 c2_max/c2_set/Weight_c2 及其 for 循环）
%   4) 雷达&目标参数整块：Num_target, delay_max, alpha_max, xi_v, Q, doppler_max*,
%      delta_tau/delta_fd/delta_range*/delta_velocity*, Max_*, velocity_*_delta_f,
%      SNR_input*, SNR_Pilot_dB, PSWR*, SNR_weak*, P_noise —— 计算后均未被使用
%   5) 仅用于物理量换算、未参与计算的整条链：c, fc, lamda, fs_up, fs, B, delta_f,
%      T, Ncp, Tcp, T_AFDM
%   6) 源头部的功能描述与脚本实际功能不符，已按实际功能重写。

clear;
close all;

%% 脉冲成形参数
L           = 2;      % 过采样倍数
beta        = 0.35;   % RRC 滚降因子
flag_pulse  = 1;      % 1: RRC 脉冲成形；0: 矩形脉冲
N_preSuffix = 20;     % 截断窗口：峰值左右各保留 N_preSuffix*L 点

%% AFDM 波形与网格参数
Nc_0                 = 128;       % 基准子载波个数
Nc                   = Nc_0;      % 子载波个数
N                    = Nc;        % 符号数
c1                   = 8/(2*Nc);  % AFDM 参数，2*N*c1 = 8
mu_4                 = 1.32;      % 16-QAM 四阶矩
L_delay_noDelayDop   = 2;         % 时延维过采样倍数
L_Doppler_noDelayDop = 2;         % 多普勒维过采样倍数

%% 生成根升余弦响应，截断后循环移位使峰值对齐
if flag_pulse == 1
    p_rrc = rcosdesign(beta, Nc, L, 'sqrt').';
    p_rrc = p_rrc(1:end-1);       % 长度 Nc*L+1 -> Nc*L
else
    p_rrc = [1; zeros(Nc*L-1, 1)];
end

idx = Nc*L/2 - N_preSuffix*L + (1 : 2*N_preSuffix*L+1);
p_rrc_trunc = zeros(Nc*L, 1);
p_rrc_trunc(idx) = p_rrc(idx);
p_rrc = circshift(p_rrc_trunc, -Nc*L/2);

%% 时延 / 多普勒网格
delay_nor_set   = 0 : 1/L_delay_noDelayDop : N*L - 1/L_delay_noDelayDop;
Doppler_nor_set = -N*L/2 : 1/L_Doppler_noDelayDop : N*L/2 - 1/L_Doppler_noDelayDop;

%% 计算理论模糊函数（tau0 = 0, v0 = 0）
% AmbigFunc_pulsedAFDM_noDelayDop_v6 = f_theoretical_AmbigFunc_pulseShaped_AFDM_v6_clean( ...
%     p_rrc, N, L, c1, mu_4, delay_nor_set, Doppler_nor_set, 0, 0);
AmbigFunc_pulsedAFDM_noDelayDop_v7 = f_theoretical_AmbigFunc_pulseShaped_AFDM_v7( ...
    p_rrc, N, L, c1, mu_4, delay_nor_set, Doppler_nor_set, 0, 0);
%% 保存（源脚本此处为注释状态；取消注释即可保存，目录不存在会自动创建）
if exist('result', 'dir') ~= 7, mkdir('result'); end
% save('./result/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_v6.mat', ...
%     'AmbigFunc_pulsedAFDM_noDelayDop_v6');
save('./result/AmbiguFunc/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_v7.mat', ...
    'AmbigFunc_pulsedAFDM_noDelayDop_v7');
