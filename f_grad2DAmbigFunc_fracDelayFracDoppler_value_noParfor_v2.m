function E_mat_total = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v2( ...
    X, P, norDelay_strong, norDopper_strong, para, echo_rows, echo_cols, cfg_trial)
%% 版本信息
% 2025/8/3 V1.0
% 2026/7/20 V2.0：保持原二维非整数循环移位定义，将逐列、逐行循环改为
%                  沿指定维度的矩阵 FFT 运算，并规范变量命名和注释。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2025  BUAA
%==========================================================================
% Designer：Ni Yuanhan
% Design Date：202500803
% Modifier：Ni Yuanhan
% Modifying Date：202500803
% Function：f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v2
% Description：根据推导的模糊函数解析式，利用二维循环移位，生成非整数
%              时延和非整数多普勒对应的目标回波采样值。
% Programming Idea：分别沿时延维和多普勒维进入频域，在中心化频率 bin 上
%                   乘线性相位，再通过逆 FFT 得到二维非整数循环移位结果。
% Subordination：无
%==========================================================================
% Modification：
%   1. 保留 V1 的函数接口、相位符号、fftshift 顺序和实部输出方式。
%   2. 取消逐列和逐行循环，直接沿矩阵第 1、2 维完成相同计算。
%   3. 保留原函数末尾被注释的搜索区域截取代码，便于需要时恢复。
%
% 输入参数：
%   X                         - 提前离线产生的无时延、无多普勒模糊函数矩阵
%   P                         - 功率
%   N                         - 子载波数量
%   L                         - 时域过采样倍数
%   L_delay_noDelayDop        - 理论模糊函数时延域过采样倍数
%   L_Doppler_noDelayDop      - 理论模糊函数多普勒域过采样倍数
%   norDelay_strong           - 目标回波对应的归一化时延
%   norDopper_strong          - 目标回波对应的归一化多普勒
%
% 输出参数：
%   E_mat_total               - 完成非整数时延和多普勒循环移位后的实值二维响应
%
% Others：
%   N、L、L_delay_noDelayDop 和 L_Doppler_noDelayDop 在当前完整矩阵输出路径
%   中不直接参与运算；保留这些参数是为了兼容 V1 接口及末尾的搜索区域截取。
% Record：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% 1. 施加目标功率因子
    % 对于上采样的脉冲成形，为保证脉冲成形前后信号功率一致，信号幅度需
    % 乘 sqrt(L)，最终造成平方模糊函数包含 L^2 因子。这里沿用 V1 的输入
    % 定义，直接将目标功率因子乘到离线理论模糊函数矩阵上。
    AmbigFunc_noDelayDoppler_scaled = P * X;

    [num_delay, num_doppler] = size(AmbigFunc_noDelayDoppler_scaled);
    AmbigFunc_fracDelayDoppler_acc = zeros([length(echo_rows), length(echo_cols)]);
    for sym_index = 1:para.Nsym
        %% 2. 沿时延维实现非整数循环移位
        % V1 对每一列执行：
        %   ifft(fftshift(fftshift(fft(column)) .* delay_phase))
        % 此处显式指定第 1 维，对所有列同时执行完全相同的运算。
        delay_bins_centered = ((0:num_delay - 1) - num_delay / 2).';
        delay_phase = exp(-1i * 2 * pi / num_delay * norDelay_strong(sym_index) * delay_bins_centered);

        delay_spectrum_centered = fftshift(fft(AmbigFunc_noDelayDoppler_scaled, [], 1), 1);
        delay_spectrum_shifted = fftshift(delay_spectrum_centered .* delay_phase, 1);
        AmbigFunc_fracDelay = ifft(delay_spectrum_shifted, [], 1);

        %% 3. 沿多普勒维实现非整数循环移位
        % V1 对时延移位结果的每一行执行：
        %   ifft(fftshift(fftshift(fft(row)) .* doppler_phase))
        % 此处显式指定第 2 维，对所有行同时执行，避免逐行循环。
        doppler_bins_centered = (0:num_doppler - 1) - num_doppler / 2;
        doppler_phase = exp(-1i * 2 * pi / num_doppler * norDopper_strong * doppler_bins_centered);

        doppler_spectrum_centered = fftshift(fft(AmbigFunc_fracDelay, [], 2), 2);
        doppler_spectrum_shifted = fftshift(doppler_spectrum_centered .* doppler_phase, 2);
        AmbigFunc_fracDelayDoppler = ifft(doppler_spectrum_shifted, [], 2);

        % 根据前期补偿的徙动做循环移位,先降采样 -> 再循环移位
        AmbigFunc_fracDelayDoppler = AmbigFunc_fracDelayDoppler(echo_rows, echo_cols);
        delay_shift = -cfg_trial.migration_compensate_bin(1, sym_index);
        doppler_shift = -cfg_trial.migration_compensate_bin(2, sym_index);
        AmbigFunc_fracDelayDoppler_cir = circshift(AmbigFunc_fracDelayDoppler, [delay_shift, doppler_shift]);

        AmbigFunc_fracDelayDoppler_acc = AmbigFunc_fracDelayDoppler_acc + AmbigFunc_fracDelayDoppler_cir;
    end
    AmbigFunc_fracDelayDoppler_acc = AmbigFunc_fracDelayDoppler_acc / para.Nsym;
    %% 4. 输出实值二维目标响应
    % 理论平方模糊函数为实值；FFT 数值计算可能留下机器精度量级的虚部。

    E_mat_total = real(AmbigFunc_fracDelayDoppler_acc);

end
