function [Tx_signal, Tx_signal_ref, qam_frame, data_bits, qam_symbols] = f_AFDM_TxSignal_v2(para, g_ps, P_rrc)
%%  版本信息
%%
%%创建日期：2026/06/16 V1.0
%%修改日期：2026/11/04 V2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：生成随机比特、完成 QAM 调制，并生成 AFDM 发射波形和匹配滤波参考信号
% Modification:去掉没用的"mod_para"参数
% 输入参数：
%   para           - AFDM、脉冲成形和雷达参数结构体
%   g_ps           - 发射端脉冲成形滤波器
%   P_rrc          - 有效 AFDM 符号的等效脉冲成形矩阵
%
% 输出参数：
%   tx_signal                - 带 CP/保护区并经过脉冲成形后的时域发射信号
%   Tx_signal_ref            - 去 CP 有效 AFDM 符号对应的匹配滤波参考矩阵，每列对应一个 AFDM 符号
%   qam_frame                - 随机 QAM 数据帧，每列对应一个 AFDM 符号
%   data_bits                - 当前 trial 随机生成的原始比特
%   qam_symbols              - 当前 trial 的 QAM 复符号序列
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;

    %% 1. QAM 数据生成
    % 每个 Monte Carlo trial 都重新生成随机 QAM 数据。
    % 每个 QAM 复符号由 log2(para.Mod) 个比特映射得到。
    bit_num = para.Nc;
    data_bits = randi([0 1], bit_num * para.Nsym * log2(para.Mod), 1);

    % data_bits_ = randi([0 1], bit_num, 1);
    % data_bits = repmat(data_bits_, para.Nsym * log2(para.Mod), 1);
    % 比特分组后转十进制索引，再映射到单位平均功率 QAM 星座点。
    data_bits_reshape = reshape(data_bits, log2(para.Mod), []).';
    data_symbols = bi2de(data_bits_reshape);
    qam_symbols = qammod(data_symbols, para.Mod, 'UnitAveragePower', true).';

    % 按 fs 和 Tb 做重复采样，并 reshape 成逐列 IDAFT 所需的帧结构。这里做上菜样为了做扩频
    kron_vec = ones(1, 1);
    qam_upsamp = kron(qam_symbols, kron_vec);
    qam_frame = reshape(qam_upsamp, length(qam_upsamp) / para.Nsym, []);

    %% 2. AFDM 调制（IDAFT）
    % 逐列 IDAFT：c2 相位加权、归一化 IFFT、c1 相位加权。
    % 相位列向量等价于原 DAFT 对角矩阵共轭转置后的对角元素。
    chirp_c2 = conj(exp(-1i * 2 * pi * para.c2 * ((0:para.M - 1).').^2));
    chirp_c1 = conj(exp(-1i * 2 * pi * para.c1 * ((0:para.N - 1).').^2));
    transmit_data_c2 = kron(ones(1, para.Nsym), chirp_c2) .* qam_frame;
    ifft_data = sqrt(para.N) * ifft(transmit_data_c2, para.N, 1);
    transmit_data_noCP_AFDM = kron(ones(1, para.Nsym), chirp_c1) .* ifft_data;

    %% 3. 添加 CP 和前后保护区
    % 为每个 AFDM 符号添加前保护、CP 和后保护，避免脉冲成形过渡影响有效区。
    transmit_data_CP_preSuffix_AFDM = [ ...
        transmit_data_noCP_AFDM(end - para.Ncp - para.N_preSuffix + 1:end - para.Ncp, :); ...
        transmit_data_noCP_AFDM(end - para.Ncp + 1:end, :); ...
        transmit_data_noCP_AFDM; ...
        transmit_data_noCP_AFDM(1:para.N_preSuffix, :)];

    %% 4. 串行化、上采样和脉冲成形
    % 串行化、L 倍上采样并通过脉冲成形滤波器，得到实际送入时变信道的发射信号。
    x_sigCP_preSuffix = reshape(transmit_data_CP_preSuffix_AFDM, [], 1).';
    x_upSamCP_preSuffix = reshape([sqrt(para.L) * x_sigCP_preSuffix; ...
        zeros(para.L - 1, samples_per_symbol * para.Nsym)], [], 1);
    % 脉冲成形
    x_pulseShaped_CP_preSuffix = conv(g_ps, x_upSamCP_preSuffix);
    Tx_signal = x_pulseShaped_CP_preSuffix;

    %% 5. 匹配滤波参考信号生成
    % 为接收端时域匹配滤波生成无 CP 有效符号的 pulse-shaped 参考信号。
    x_upSam = kron(sqrt(para.L) * transmit_data_noCP_AFDM, [1; zeros(para.L - 1, 1)]);
    Tx_signal_ref = P_rrc * x_upSam;
end
