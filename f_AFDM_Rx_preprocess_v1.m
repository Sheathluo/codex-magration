function X_PS_remove_GPS = f_AFDM_Rx_preprocess_v1(Rx_signal_AFDM_PS_GPS, cfg, para)
%%  版本信息
%%
%%2026/06/16 V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：截取接收信号，去成形过渡带，并提取每个 AFDM 符号有效区域及CP
%
% 输入参数：
%   Rx_signal_AFDM_PS_GPS   - 时变信道输出的完整接收信号，可能包含最大时延扩展
%   cfg                     - 接收配置结构体，主要使用 cfg.index_start 指定截取起点
%   para                    - AFDM 与脉冲成形参数结构体
%   Len_sig                 - 主函数计算的接收截取长度，与发射有效观测帧长度对齐
%
% 输出参数：
%   X_PS_remove_GPS - 接收端去前后保护区后的有效 AFDM 符号矩阵
%                                    行为过采样后的有效采样点，列为 AFDM 符号索引
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;
    samples_per_symbol_L = samples_per_symbol * para.L;
    % 去掉过信道补的0
    Rx_signal_AFDM = Rx_signal_AFDM_PS_GPS(cfg.index_start + 1 : cfg.index_start + para.len_signal);

    % 去掉卷积脉冲成形带来的首尾过渡段，再按符号长度重排成帧。
    start_trunc = para.N_tapFiler * para.L / 2 + 1;
    end_trunc = para.N_tapFiler * para.L / 2 + samples_per_symbol_L * para.Nsym;
    Rx_signal_AFDM_trunc = Rx_signal_AFDM(start_trunc:end_trunc);

    rx_frame = reshape(Rx_signal_AFDM_trunc, samples_per_symbol_L, para.Nsym);
    % 保留 CP + 有效 AFDM 符号区，去除额外前后保护采样。
    valid_start = (para.Ncp + para.N_preSuffix) * para.L + 1;
    valid_idx = valid_start : valid_start + (para.N + para.Ncp) * para.L - 1;
    X_PS_remove_GPS = rx_frame(valid_idx, :);

end
