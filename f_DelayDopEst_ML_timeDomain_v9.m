function [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v9(waveform_index,factor_delay_refine, ...
    factor_Doppler_refine,Tx_signal,Rx_signal,Min_delay,Max_delay,Min_nor_Doppler,Max_nor_Doppler, ...
    delta_fd,delta_f,fs,L,Num_target,N,Ncp,Nsym)
%%  版本信息
%%
%%2024/11/14 V2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2024  BUAA
%======================================================
% Designer：Ni Yuanhan
% Design Date：20241114
% Modifier：Ni Yuanhan
% Modifying Date：20241114
% Function
% Description：在时域，利用最大似然估计算法， 估计目标时延（整数）、多普勒和增益
% Programming
% Idea：
% Input parameter:  Tx_signal - 发送波形
%                   Rx_signal - 接收波形,
%                   Max_delay - 需要搜索的时延范围[0,Max_delay]，
%                   Max_Doppler - 需要搜索的时延范围[-Max_Doppler,Max_delay],
%                   delta_fd - 多普勒分辨率
%                   Ncp - CP点数
%                   Nsym - 列数
%                   transmit_data - 发送数据
% Output parameter: try_catch - 错误指示标识,
%                   delay_est - 估计的时延,
%                   Doppler_est - 估计的多普勒,
%                   chan_coef_est - 估计的信道增益
% Subordination：无
% Others：
%======================================================
% Modification： v2在v1的基础上，增加了根据不同波形类型，自适应的调整去CP的方式
%                v3在v2的基础上，增加了Doppler搜索因子可选的功能
%                v6在v5的基础上，增加了利用循环移位进行非整数时延搜索的功能
%                v9在v6的基础上，更改了匹配滤波逻辑，将时延循环移位从作用在参考信号上转变为作用在接收信号上
% Record：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% *********** 粗搜索：先搜索整数时延-多普勒分辨单元 **************%
% Nor_delay_roughSearch_set = Min_delay:Max_delay;
% Max_Doppler = Max_nor_Doppler * delta_f;
% Num_Doppler_res = ceil(Max_Doppler /delta_fd);
% Nor_Doppler_roughSearch_set = (-Num_Doppler_res+1:1:Num_Doppler_res) * delta_fd / delta_f ;
%
% Len_delay = length(Nor_delay_roughSearch_set);
% Len_Doppler = length(Nor_Doppler_roughSearch_set);
%
% ML_value = zeros(Len_delay,Len_Doppler);
%
% % 计算统计量
% for index_delay = 1:Len_delay
%     nor_delay = Nor_delay_roughSearch_set(index_delay);
%     for index_Doppler = 1:Len_Doppler
%         nor_Doppler = Nor_Doppler_roughSearch_set(index_Doppler);
%
%         % 本地生成回波
%         chan_coef_est = 1;
%         [try_catch, Rx_signal_copy_cp] = f_passSensingChannel_v1(Tx_signal,N,Ncp,Num_target,chan_coef_est, nor_delay, nor_Doppler, fs);
%         % 去CP
%         if strcmp(waveform_index, 'OFDM') || strcmp(waveform_index, 'AFDM') % 全CP OFDM，二维块
%             Rx_signal_copy_cp_2D = reshape(Rx_signal_copy_cp,N+Ncp,Nsym);
%             Rx_signal_copy_2D = Rx_signal_copy_cp_2D(Ncp+1:end,:);
%             Rx_signal_copy = reshape(Rx_signal_copy_2D,[],1);
%         else % RCP AFDM，一维向量
%             Rx_signal_copy = Rx_signal_copy_cp(Ncp+1:end);
%         end
%
%         % 计算统计量
%         ML_value(index_delay,index_Doppler) = abs(Rx_signal_copy'*Rx_signal).^2 / (norm(Rx_signal_copy).^2);
%
%     end
% end
%
% % 搜索峰值
% [value_row_max, I_row_set_max] = max(ML_value);
% [value_max, I_col_max] = max(value_row_max);
% I_row_max = I_row_set_max(I_col_max);
% % ---------------- 结束：粗搜索过程 -------------------------------------------%
%
% %% ************** 小范围精搜索 **************************************%

% % factor_Doppler_refine = 0.001;
% del_norDoppler = factor_Doppler_refine*delta_fd / delta_f;
%
% I_delay_begin = I_row_max-1;
% I_delay_end = I_row_max+1;
% I_Doppler_begin = I_col_max-1;
% I_Doppler_end = I_col_max+1;
%
% % 避免索引越界
% if I_delay_begin < 1
%     I_delay_begin = 1;
% end
% if I_delay_end > Len_delay
%     I_delay_end = Len_delay;
% end
% if I_Doppler_begin < 1
%     I_Doppler_begin = 1;
% end
% if I_Doppler_end > Len_Doppler
%     I_Doppler_end = Len_Doppler;
% end

del_norDelay = factor_delay_refine;
del_norDoppler = factor_Doppler_refine;

Nor_delay_refineSearch_set = Min_delay:del_norDelay:Max_delay;
Nor_Doppler_refineSearch_set = Min_nor_Doppler:del_norDoppler:Max_nor_Doppler;

Len_delay_refine = length(Nor_delay_refineSearch_set);
Len_Doppler_refine = length(Nor_Doppler_refineSearch_set);
ML_value_refine = zeros(Len_delay_refine,Len_Doppler_refine);
% 计算统计量
for index_delay = 1:Len_delay_refine
    nor_delay_refine = Nor_delay_refineSearch_set(index_delay);
    for index_Doppler = 1:Len_Doppler_refine
        nor_Doppler_refine = Nor_Doppler_refineSearch_set(index_Doppler);

        % 本地生成回波
        chan_coef_est = 1;
        % 先尝试时延副本
        % v6
%         Tx_signal_copy_fracDelay = ifft(fftshift(fftshift(fft(Tx_signal)).*exp(-1i*2*pi./(length(Tx_signal)).*(nor_delay_refine)*((0:length(Tx_signal)-1)-length(Tx_signal)/2).')));
        % v9
        Tx_signal_copy_fracDela = Tx_signal;
        %         figure
%         plot(real(Tx_signal),'-r*'); hold on;
%         plot(real(Tx_signal_copy_fracDelay),'--bo'); hold on;
%         grid on;
%         legend('原信号','循环移位后信号');

        Rx_signal_copy = conj(Tx_signal) .* exp(-1i.*2.*pi.*nor_Doppler_refine/(N*L).*(0:N*L-1).');
%         [try_catch, Rx_signal_copy] = f_passSensingChannel_v2(Tx_signal,N,Ncp,Num_target,chan_coef_est, nor_delay_refine, nor_Doppler_refine, fs, L);
        % 去CP
%         if strcmp(waveform_index, 'OFDM') || strcmp(waveform_index, 'AFDM') % 全CP OFDM，二维块
%             Rx_signal_copy_cp_2D = reshape(Rx_signal_copy_cp,N+Ncp,Nsym);
%             Rx_signal_copy_2D = Rx_signal_copy_cp_2D(Ncp+1:end,:);
%             Rx_signal_copy = reshape(Rx_signal_copy_2D,[],1);
%         else % RCP AFDM，一维向量
%             Rx_signal_copy = Rx_signal_copy_cp(Ncp+1:end);
%         end

%         Rx_signal_fracDelay = ifft(fftshift(fftshift(fft(Rx_signal)).*exp(-1i*2*pi./(length(Rx_signal)).*(-nor_delay_refine)*((0:length(Rx_signal)-1)-length(Rx_signal)/2).')));
        Rx_signal_fracDelay = Rx_signal(1+nor_delay_refine:N*L+nor_delay_refine,1);

        % 计算统计量
%         ML_value_refine(index_delay,index_Doppler) = abs(Rx_signal_copy'*Rx_signal).^2 ;
        ML_value_refine(index_delay,index_Doppler) = 1/(N*L) * (Rx_signal_copy.'*Rx_signal_fracDelay);
%         ML_value_refine(index_delay,index_Doppler) = Rx_signal'*Rx_signal_copy ;



    end
end

% % 搜索峰值
% [value_row_max_refine, I_row_set_max_refine] = max(ML_value_refine,[],1);
% [value_max_refine, I_col_max_refine] = max(value_row_max_refine);
% I_row_max_refine = I_row_set_max_refine(I_col_max_refine);
%
% % 计算时延和多普勒估计值
% nor_delay_est = Nor_delay_refineSearch_set(I_row_max_refine);
% nor_Doppler_est = Nor_Doppler_refineSearch_set(I_col_max_refine);
%
% delay_est = nor_delay_est/fs;
% Doppler_est = nor_Doppler_est*delta_f;

% --------------------------------------------------%

end