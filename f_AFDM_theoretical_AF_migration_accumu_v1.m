function AFDM_AF_migration_accumulation = f_AFDM_theoretical_AF_migration_accumu_v1(...
    para, tau0, v0, delay_bins, doppler_bins,delay_scale, doppler_scale, AFDM_AF_freq, gamma,...
    delay_range, doppler_range,index_delay_set_fracEst, index_Doppler_set_fracEst)

%% 版本信息
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026  BUAA
%==========================================================================
% Designer：Luo Tao
% Design Date：20260718
% Function：f_AFDM_theoretical_AF_migration_accumu_v1
% Description：根据粗估计的多普勒和时延生成理论模糊函数
%==========================================================================
% Modification：
%
% 输入参数：
%   para                    接收回波信号的符号数量
%   tau0                    当前初始归一化时延估计值
%   v0                      当前归一化多普勒估计值
%   fc                      载波频率
%   L                       过采样倍数
%   delay_bins              时延网格
%   doppler_bins            多普勒网格
% 输出参数：
%   response_migration      有徙动的理论模糊函数的叠加
%
% Record：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    AFDM_sym_index = 0:para.Nsym-1;
    % 根据估计出的时延tau0和多普勒v0构造理论模糊函数
    tau_migration = tau0 + AFDM_sym_index.' .* gamma * v0;
    response_migration = zeros(size(AFDM_AF_freq));
    doppler_phase = exp(doppler_scale * v0 * doppler_bins);
    for sym_index = 1:para.Nsym
        delay_phase = exp(delay_scale * tau_migration(sym_index) * delay_bins).';
        response_sym = real(ifft(fft(ifft(delay_phase .* AFDM_AF_freq, [], 1), [], 2) .* doppler_phase, [], 2));
        response_migration = response_migration + response_sym;
    end
    response_migration = response_migration / para.Nsym;
    response_migration = response_migration(delay_range, doppler_range);
    AFDM_AF_migration_accumulation = response_migration(index_delay_set_fracEst, index_Doppler_set_fracEst);

end
