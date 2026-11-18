function E_mat_total = f_grad2DAmbigFunc_fracDelayFracDoppler_value_v2(X,N,L,L_delay_noDelayDop,L_Doppler_noDelayDop, ...
    P_sig_strong,norDelay_strong,norDopper_strong)
%%  版本信息
%%
%%2025/8/3 V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2025  BUAA
%======================================================
% Designer：Ni Yuanhan
% Design Date：202500803
% Modifier：Ni Yuanhan
% Modifying Date：202500803
% Function
% Description：根据推导的模糊函数解析式，利用二维循环移位，生成非整数时延和非整数多普勒对应的目标回波采样值
% Programming
% Idea：
% Input parameter:  X - 提前离线产生的无时延、无多普勒的模糊函数
%                   N - 子载波数量
%                   L - 时域过采样倍数
%                   L_delay_noDelayDop - 时延域过采样倍数
%                   L_Doppler_noDelayDop - 多普勒域过采样倍数
%                   max_norDelay - 归一化时延（与过采样后的采样率归一化）的范围
%                   max_norDoppler - 归一化多普勒（与子载波归一化）的范围
%                   P_sig_strong - 目标回波对应的功率因子
%                   norDelay_strong - 目标回波对应的归一化时延（与过采样后的采样率归一化）
%                   norDopper_strong - 目标回波对应的归一化多普勒（与子载波归一化）
% Subordination：无
% Others：
%======================================================
% Modification：
% Record：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 对于上采样的脉冲成形，为了保证，脉冲成形前后信号功率一致，信号幅度需乘sqrt(L)，最终造成平方模糊函数成L^2
    AmbigFunc_pulsedAFDM_noDelayDop_pulsed = P_sig_strong * X;
    % 尝试非整数循环移位（利用FFT实现）
    AmbigFunc_pulsedAFDM_fracDelay = zeros(size(AmbigFunc_pulsedAFDM_noDelayDop_pulsed,1),size(AmbigFunc_pulsedAFDM_noDelayDop_pulsed,2));
    AmbigFunc_pulsedAFDM_fracDelayDoppler = zeros(size(AmbigFunc_pulsedAFDM_noDelayDop_pulsed,1),size(AmbigFunc_pulsedAFDM_noDelayDop_pulsed,2));
    for m = 1:size(AmbigFunc_pulsedAFDM_noDelayDop_pulsed,2)
        E_mat_strong_noDelayDop_col = AmbigFunc_pulsedAFDM_noDelayDop_pulsed(:,m);

        AmbigFunc_pulsedAFDM_fracDelay_tmp = fftshift(fftshift(fft(E_mat_strong_noDelayDop_col)).*exp(-1i*2*pi./(length(E_mat_strong_noDelayDop_col)).*norDelay_strong*((0:length(E_mat_strong_noDelayDop_col)-1)-length(E_mat_strong_noDelayDop_col)/2).'));
        AmbigFunc_pulsedAFDM_fracDelay(:,m) = ifft( AmbigFunc_pulsedAFDM_fracDelay_tmp);
    end
    for k = 1:size(AmbigFunc_pulsedAFDM_fracDelay,1)
        E_mat_strong_intDop_row = AmbigFunc_pulsedAFDM_fracDelay(k,:);
    %             E_mat_strong_intFracDop_tmp = fft(E_mat_strong_intDop_row).*exp(-1i*2*pi./(length(E_mat_strong_intDop_row)).*frac_norDopper_strong.*L_Doppler_noDelayDop.*(0:length(E_mat_strong_intDop_row)-1));
        E_mat_strong_intFracDop_tmp = fftshift(fftshift(fft(E_mat_strong_intDop_row)).*exp(-1i*2*pi./(length(E_mat_strong_intDop_row)).*norDopper_strong*((0:length(E_mat_strong_intDop_row)-1)-length(E_mat_strong_intDop_row)/2)));
        AmbigFunc_pulsedAFDM_fracDelayDoppler(k,:) = ifft( E_mat_strong_intFracDop_tmp);
    end
    E_mat_total = real(AmbigFunc_pulsedAFDM_fracDelayDoppler);
%     E_mat_total = real(AmbigFunc_pulsedAFDM_fracDelayDoppler(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2));
end