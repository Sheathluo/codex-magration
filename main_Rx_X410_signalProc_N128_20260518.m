%%  版本信息
%%
%% 2025/04/01 V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2025  BUAA
%==========================
% Designer：Ni Yuanhan
% Design Date：20250626
% Modifier：Ni Yuanhan
% Modifying Date：20250805
% Function Description：对比OFDM和AFDM在强、弱目标条件下的测速误差。
% Modification：v11:在v10的基础上，优化代码，适配X410等软件无线电实验
%               v12:在v11的基础上，增大采样率，适配X410等软件无线电实验
%               20260518在20260511基础上增加了先粗补偿距离徙动
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
rng(1); % 固定随机数种子，保障算法排查时结果绝对可复现
tic

%% =========== 1. 参数加载与初始化 ===========
para_path = './data/para_set_pulsedAFDM_N128_ISAC_time_20251124_X410.mat';
if ~exist(para_path, 'file')
    error('未找到参数文件: %s，请先运行生成脚本。', para_path);
end
load(para_path, 'para_set');

L           = para_set.L;
beta        = para_set.beta;
L_delay     = para_set.L_delay;             % 原始加载的 L_delay
L_Doppler   = para_set.L_Doppler;
flag_pulse  = para_set.flag_pulse;
flag_RRC    = para_set.flag_RRC;
N_preSuffix = para_set.N_preSuffix;
N_tapFiler  = para_set.N_tapFiler;
fs_up       = para_set.fs_up;
fs          = para_set.fs;
c           = para_set.c;
fc          = para_set.fc;
Nc          = para_set.Nc;
Nsym        = para_set.Nsym;
Mod         = para_set.Mod;
B           = para_set.B;
Ncp         = para_set.Ncp;
c1          = para_set.c1;
c2          = para_set.c2;

c = 3e8;
L_delay = 1; % 覆盖加载的参数
factor_delay_refine = 1/L_delay;
factor_Doppler_refine = 1/L_Doppler;

% 部分空间搜索参数
max_norDelay = Ncp - 1;
max_norDoppler = 15;
Min_delay = 0;
Max_delay = max_norDelay * L - 1;
Min_nor_Doppler = -max_norDoppler * L / 2;
Max_nor_Doppler = max_norDoppler * L / 2 - 1/L_Doppler;

ML_value_refine_set = zeros((Max_delay-Min_delay)+1, (Max_nor_Doppler-Min_nor_Doppler)/factor_Doppler_refine+1, Nsym);
lamda = c/fc; % 载波波长
delta_f = B/Nc;
T = 1/delta_f; % s,码元持续时间
Tcp = Ncp/B;   % s,循环前缀时间
T_AFDM = T+Tcp; % s,整个AFDM符号持续时间
inter_time = Nsym*T_AFDM;

N_delay = Nc;
N_Doppler = Nc;
N = Nc;
M = Nc;

%% =========== 2. 雷达 & 目标参数 ===========
Num_target = 2; % 目标（多径）数目
delay_max = Num_target;
alpha_max = 2; % 被delta_f归一化的多普勒频率
xi_v = 4; % 等效信道宽度扩展
Q = (delay_max + 1)*(2*(alpha_max + xi_v)+1) - 1;

delta_tau = T/Nc; % 时延分辨率
delta_fd = 1/(T_AFDM*Nsym); % 多普勒分辨率
delta_range_com = delta_tau*c; % 距离分辨率（通信，单程）
delta_velocity_com = delta_fd/fc*c; % 速度分辨率（通信，单程）m/s
delta_range_rad = delta_tau*c/2; % 距离分辨率（雷达，双程）
delta_velocity_rad = delta_fd/fc*c/2; % 速度分辨率（雷达，双程）m/s

Max_unamb_range = c*T/2;
Max_ISI_free_range = c*Tcp/2;
Max_unamb_velocity = c/(2*fc*T_AFDM);
Max_ICI_free_velocity = B*c/(20*Nc*fc);

%% =========== 3. 读取外部文件与理论数据 ===========
% 读取理论模糊函数
ambig_path = './result/AmbiguFunc/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v6_20251022.mat';
if exist(ambig_path, 'file')
    load(ambig_path, 'AmbigFunc_pulsedAFDM_noDelayDop');

else
    % error('未找到理论模糊函数文件: %s', ambig_path);
    % 生成根升余弦滤波器响应
    if flag_pulse == 1
        p_rrc = rcosdesign(beta,Nc,L,'sqrt').';
        p_rrc = p_rrc(1:end-1);
    else
        p_rrc = [1; zeros(Nc*L-1,1)];
    end

    % 峰值左右保留N_preSuffix*L个数
    p_rrc_trunc = zeros(Nc*L,1);
    p_rrc_trunc(Nc*L/2-N_preSuffix*L+1:Nc*L/2-N_preSuffix*L+(2*N_preSuffix*L+1)) = p_rrc(Nc*L/2-N_preSuffix*L+1:Nc*L/2-N_preSuffix*L+(2*N_preSuffix*L+1));

    % 循环移位
    p_rrc = circshift(p_rrc_trunc, -Nc*L/2);

    mu_4 = 1.32; % 16-QAM
    delay_nor_set = 0:1/L_delay:N*L-1/L_delay;
    Doppler_nor_set = -N*L/2:1/L_Doppler:N*L/2-1/L_Doppler;
    AmbigFunc_pulsedAFDM_noDelayDop = f_theoretical_AmbigFunc_pulseShaped_AFDM_v7(p_rrc, N, L, c1, ...
        mu_4, delay_nor_set, Doppler_nor_set, 0, 0);
    save ./result/AmbiguFunc/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v6_20251022.mat AmbigFunc_pulsedAFDM_noDelayDop

end
AmbigFunc_pulsedAFDM_noDelayDop = AmbigFunc_pulsedAFDM_noDelayDop / (N^2);

% 辅助函数：安全读取复数 dat 文件
read_complex_dat = @(f_re, f_im) ...
    complex(fread(fopen(f_re,'r'), 'float32') * (fclose(fopen(f_re,'r'))==0), ...
            fread(fopen(f_im,'r'), 'float32') * (fclose(fopen(f_im,'r'))==0));

% 读取原始发送信号
Tx_re_path = './data/Tx_pulsedAFDM_N128_ISAC_time_20251124_X410_re.dat';
Tx_im_path = './data/Tx_pulsedAFDM_N128_ISAC_time_20251124_X410_im.dat';
Tx_signal_tmp = read_complex_dat(Tx_re_path, Tx_im_path);
Tx_signal = Tx_signal_tmp;

% 读取参考信号
Ref_re_path = './data/Tx_ConfencedSignal_pulsedAFDM_N128_ISAC_time_20251124_X410_re.dat';
Ref_im_path = './data/Tx_ConfencedSignal_pulsedAFDM_N128_ISAC_time_20251124_X410_im.dat';
Ref_signal_tmp = read_complex_dat(Ref_re_path, Ref_im_path);
Tx_signal_pulsedAFDM_eff = reshape(Ref_signal_tmp(1:Nc*L*Nsym), Nc*L, Nsym);

% 读取接收信号
Rx_re_path = './result/Rx_simSensing_pulsedAFDM_N128_ISAC_time_20251124_X410_re.dat';
Rx_im_path = './result/Rx_simSensing_pulsedAFDM_N128_ISAC_time_20251124_X410_im.dat';
Rx_signal_AFDM_tmp = read_complex_dat(Rx_re_path, Rx_im_path);

% figure;
% plot(abs(Rx_signal_AFDM_tmp));
% title('原始接收信号包络');

%% =========== 4. 截取有效数据与预处理 ===========
index_start = 0;
Len_sig = (Nc+Ncp+2*N_preSuffix)*L*Nsym + N_preSuffix*L;

% 【边界保护】：防止由于截取长度过大导致索引越界崩溃
if index_start + Len_sig > length(Rx_signal_AFDM_tmp)
    error('请求截取的数据长度超出了实际接收到的数据总长！请检查 Nsym 或 index_start。');
end
Rx_signal_AFDM = Rx_signal_AFDM_tmp(index_start+1 : index_start+Len_sig);

Len_tx = length(Tx_signal);
% f_tx = fs_up/Len_tx*(0:Len_tx-1);
Len_rx = length(Rx_signal_AFDM);
F_Rx_signal_AFDM = abs(fftshift(fft(Rx_signal_AFDM)));
F_Tx_signal = abs(fftshift(fft(Tx_signal)));
% f_rx = fs_up/Len_rx*(0:Len_rx-1);
f_rx = (-Len_tx/2:(Len_tx/2-1))*(fs/Len_tx);

% figure
% plot(f_rx, 20*log10(F_Rx_signal_AFDM/max(F_Rx_signal_AFDM)), '--r'); hold on;
% plot(f_rx, 20*log10(F_Tx_signal/max(F_Tx_signal)), '-bo');
% ylim([-150,5]);
% legend('Rx', 'Tx');
% title('收发信号频谱对比');

% 去掉滤波器前后过渡带
Rx_signal_AFDM_trunc = Rx_signal_AFDM(N_tapFiler*L/2+1 : N_tapFiler*L/2+(Nc+Ncp+2*N_preSuffix)*Nsym*L);
% 移出每个符号的保护间隔和CP
X_pulseShaped_remove_preSuffix = zeros((N+Ncp)*L, Nsym);
for n_sym = 1:Nsym
    idx_start = (Nc+Ncp+2*N_preSuffix)*L*(n_sym-1) + (Ncp+N_preSuffix)*L + 1;
    X_pulseShaped_remove_preSuffix(:,n_sym) = Rx_signal_AFDM_trunc(idx_start : idx_start+(Nc+Ncp)*L-1);
end

% X_pulseShaped_CP_preSuffix_rep = reshape(Rx_signal_AFDM_trunc, (Nc+Ncp+2*N_preSuffix)*L, Nsym);
% X_pulseShaped_remove_CP_preSuffix = X_pulseShaped_CP_preSuffix_rep((Ncp+N_preSuffix)*L+1 : (Ncp+N_preSuffix)*L+Nc*L, :);

% 计算噪声与 SNR
P_signal_noise_cal = mean(abs(Rx_signal_AFDM_trunc).^2);
index_start_noise = 1e3;
if index_start_noise + Len_sig <= length(Rx_signal_AFDM_tmp)
    noise_received = Rx_signal_AFDM_tmp(index_start_noise+Len_sig : end);
    P_noise_cal = mean(abs(noise_received).^2);
    SNR_dB_cal = 10*log10((P_signal_noise_cal-P_noise_cal)/P_noise_cal);
else
    P_noise_cal = mean(abs(Rx_signal_AFDM_tmp(1:100)).^2); % 备用噪声计算
end

%% =========== 5. 接收端时域匹配滤波 ===========
Num_target_est = 1;
waveform_index = 'AFDM';

Nor_Doppler_set = Min_nor_Doppler : factor_Doppler_refine : Max_nor_Doppler;
% Nor_Doppler_set = [12];
est_doppler_len = length(Nor_Doppler_set);
ML_value_refine_set_vactor = zeros([size(ML_value_refine_set,[1,2]), est_doppler_len]);
X_pulseShaped_remove_preSuffix_circ = zeros(size(X_pulseShaped_remove_preSuffix));
% disp('正在执行时域匹配滤波...');
disp('正在执行粗估计并时域匹配滤波...');
% 粗估计并补偿距离徙动
for dopple_index = 1:est_doppler_len
    velocity = Nor_Doppler_set(dopple_index) .* delta_f * c / (2*fc);
    migration_delta_tau_norm = (1:Nsym) * 2 .* (velocity'.*L*(Nc+Ncp+2*N_preSuffix)/c);   % 用估计的速度计算徙动

    migration_delta_tau_time = migration_delta_tau_norm/fs;
    % 根据多普勒的正负来判断移位
    if Nor_Doppler_set(dopple_index) >= 0
        range_migration_bin = floor(migration_delta_tau_time/delta_tau);
    else
        range_migration_bin = ceil(migration_delta_tau_time/delta_tau);
    end
    for n_sym = 1:Nsym
        X_pulseShaped_remove_preSuffix_circ(:, n_sym) = circshift(X_pulseShaped_remove_preSuffix(:, n_sym), ...
            -range_migration_bin(n_sym));
    end
    % 时域匹配滤波
    parfor n_sym = 1:Nsym
    % for n_sym = 1:Nsym
        % % v9用for循环很慢
        % [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v9(waveform_index, ...
        %     factor_delay_refine, factor_Doppler_refine, Tx_signal_pulsedAFDM_eff(:,n_sym), ...
        %     X_pulseShaped_remove_preSuffix(:,n_sym), Min_delay, Max_delay, Min_nor_Doppler, ...
        %     Max_nor_Doppler, delta_fd, delta_f, fs, L, Num_target_est, N, Ncp, Nsym);
        % % v10去掉for循环，实现加速
        [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v10(factor_delay_refine, factor_Doppler_refine, ...
            Tx_signal_pulsedAFDM_eff(:,n_sym),X_pulseShaped_remove_preSuffix_circ(:,n_sym), ...
            Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler, L, N);
        % % 补偿距离徙动
        % [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v11(factor_delay_refine, factor_Doppler_refine, ...
        %     Tx_signal_pulsedAFDM_eff(:,n_sym),X_pulseShaped_remove_preSuffix(:,n_sym), ...
        %     Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler, L, N);
        ML_value_refine_set(:,:,n_sym) = ML_value_refine;

    end
    % 非相干累积，估计一个多普勒，做一次非相干累积
    ML_value_refine_set_vactor(:,:,dopple_index) = sum(abs(ML_value_refine_set).^2, 3) ./ Nsym;
end
% % 非相干累积
% ML_value_refine_set_aver = sum(abs(ML_value_refine_set).^2, 3) ./ Nsym;
[max_val, max_idx] = max(ML_value_refine_set_vactor(:));
[~, ~, max_page] = ind2sub(size(ML_value_refine_set_vactor), max_idx);
ML_value_refine_set_aver = ML_value_refine_set_vactor(:,:,max_page);
% 相干累积：先对复数求和，再取模平方
% ML_value_refine_set_aver_co = abs(sum(ML_value_refine_set, 3)).^2 ./ (Nsym^2);

del_norDelay = factor_delay_refine;
del_norDoppler = factor_Doppler_refine;
Nor_delay_refineSearch_set = Min_delay:del_norDelay:Max_delay;
Nor_Doppler_refineSearch_set = Min_nor_Doppler:del_norDoppler:Max_nor_Doppler;

[X, Y] = meshgrid(Nor_Doppler_refineSearch_set, Nor_delay_refineSearch_set);

% 非相干累积
figure
mesh(X,Y, 10*log10(ML_value_refine_set_aver/max(max(ML_value_refine_set_aver))));
ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
zlabel('Normalized average squared MF');
set(gca, 'FontSize', 10);
axis([-inf inf -inf inf -40 2]);
title('单目标-粗补偿徙动后-实际匹配滤波结果-非相干累积');

% 限制动态范围，让图美观
% floor_limit = -40;
% db_data = 10*log10(ML_value_refine_set_aver_co/max(max(ML_value_refine_set_aver_co)));
% db_data(db_data < floor_limit) = floor_limit;
% 相干累积
% figure
% mesh(X,Y, db_data);
% ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
% xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
% zlabel('Normalized average squared MF');
% set(gca, 'FontSize', 10);
% axis([-inf inf -inf inf floor_limit-5 2]);
% title('单目标-粗补偿徙动后-实际匹配滤波结果-相干累积');

%%%%%%%%%%%%%%%%%%%%%%%不补偿徙动%%%%%%%%%%%%%%%%%%%%%%%%%
parfor n_sym = 1:Nsym
% for n_sym = 1:Nsym
    % % v9用for循环很慢
    % [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v9(waveform_index, ...
    %     factor_delay_refine, factor_Doppler_refine, Tx_signal_pulsedAFDM_eff(:,n_sym), ...
    %     X_pulseShaped_remove_preSuffix(:,n_sym), Min_delay, Max_delay, Min_nor_Doppler, ...
    %     Max_nor_Doppler, delta_fd, delta_f, fs, L, Num_target_est, N, Ncp, Nsym);
    % % v10去掉for循环，实现加速
    [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v10(factor_delay_refine, factor_Doppler_refine, ...
        Tx_signal_pulsedAFDM_eff(:,n_sym),X_pulseShaped_remove_preSuffix(:,n_sym), ...
        Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler, L, N);
    % % 补偿距离徙动
    % [ML_value_refine] = f_DelayDopEst_ML_timeDomain_v11(factor_delay_refine, factor_Doppler_refine, ...
    %     Tx_signal_pulsedAFDM_eff(:,n_sym),X_pulseShaped_remove_preSuffix(:,n_sym), ...
    %     Min_delay, Max_delay, Min_nor_Doppler, Max_nor_Doppler, L, N);
    ML_value_refine_set(:,:,n_sym) = ML_value_refine;
end

ML_value_refine_set_aver = sum(abs(ML_value_refine_set).^2, 3) ./ Nsym;
[X, Y] = meshgrid(Nor_Doppler_refineSearch_set, Nor_delay_refineSearch_set);

% 非相干累积
figure
mesh(X,Y, 10*log10(ML_value_refine_set_aver/max(max(ML_value_refine_set_aver))));
ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
zlabel('Normalized average squared MF');
set(gca, 'FontSize', 10);
axis([-inf inf -inf inf -40 2]);
title('单目标-不补偿徙动-实际匹配滤波结果-非相干累积');
%%%%%%%%%%%%%%%%%%%%%%%不补偿徙动%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========== 6. 梯度下降法估计非整数时延和多普勒 ===========
% 先将模糊函数反折 X(-t,-v)
AmbigFunc_pulsedAFDM_noDelayDop_flip = zeros(size(AmbigFunc_pulsedAFDM_noDelayDop));
for i = 0:size(AmbigFunc_pulsedAFDM_noDelayDop,1)-1
    i_flip = mod(-i, size(AmbigFunc_pulsedAFDM_noDelayDop,1));
    for k = 0:size(AmbigFunc_pulsedAFDM_noDelayDop,2)-1
        k_flip = mod(-k, size(AmbigFunc_pulsedAFDM_noDelayDop,2));
        AmbigFunc_pulsedAFDM_noDelayDop_flip(i+1, k+1) = AmbigFunc_pulsedAFDM_noDelayDop(i_flip+1, k_flip+1);
    end
end
% 反折模糊函数 X(-t,-v),优化版
% [M, N] = size(AmbigFunc_pulsedAFDM_noDelayDop);
% AmbigFunc_pulsedAFDM_noDelayDop_flip_ = AmbigFunc_pulsedAFDM_noDelayDop([1, M:-1:2], [1, N:-1:2]);

N_delay_peak = 4;
N_Doppler_peak = 4;
Num_iter_single = 1;
Num_iter_cyclic = 30;
Pfa_target = 1e-3;
P_noise_MF = P_noise_cal ./ N;
Th_noTarget = 1/Nsym * gaminv(1-Pfa_target, Nsym, P_noise_MF);
Max_num_target = Num_target;

L_delay_noDelayDop = 2;
L_Doppler_noDelayDop = 2;
L_delay_aver = 1/(Nor_delay_refineSearch_set(2)-Nor_delay_refineSearch_set(1));
L_Doppler_aver = 1/(Nor_Doppler_refineSearch_set(2)-Nor_Doppler_refineSearch_set(1));

disp('利用理论模糊函数进行参数精细估计...');
[set_nor_parameters_esti] = f_paraEstibyGradDesc_newtonOMP_multiTarget_v1(ML_value_refine_set_aver, ...
    AmbigFunc_pulsedAFDM_noDelayDop_flip, Nor_delay_refineSearch_set, Nor_Doppler_refineSearch_set, N, L, ...
    L_delay, L_Doppler, max_norDelay, max_norDoppler, N_delay_peak, N_Doppler_peak, fs_up, delta_f, Num_iter_single, ...
    Num_iter_cyclic, Th_noTarget, Max_num_target);

nor_delay_esti = set_nor_parameters_esti(1);
nor_Doppler_esti = set_nor_parameters_esti(2);
P_sig_strong_esti = set_nor_parameters_esti(3);

%%%% 对接收回波徙动精细补偿Precise estimation
% precise_velocity = nor_Doppler_esti .* delta_f * c / fc;   % 精估计后目标真实速度
% precise_migration_delta_tau_norm = (1:Nsym) * 2 .* (precise_velocity'.*L*(Nc+Ncp+2*N_preSuffix)/c);   % 用估计的速度计算徙动
% precise_migration_delta_tau_time = precise_migration_delta_tau_norm/fs;


%% =========== 7. 干扰消除与结果呈现 ===========
nor_delay_esti_strong_upSam = nor_delay_esti * L_delay_aver * (L_delay_noDelayDop/L_delay_aver);
nor_Doppler_esti_strong_upSam = nor_Doppler_esti * L_Doppler_aver * (L_Doppler_noDelayDop/L_Doppler_aver);

% 反演强目标回波
D_func_search_tmp = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor(AmbigFunc_pulsedAFDM_noDelayDop, ...
    N, L, L_delay_noDelayDop, L_Doppler_noDelayDop, P_sig_strong_esti, nor_delay_esti_strong_upSam, ...
    nor_Doppler_esti_strong_upSam);
% 从反演出的目标回波中采样，根据实际接收回波才来
echo_strong_esti = D_func_search_tmp(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop, ...
    N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);


% 从接收回波中减去强目标回波，得到弱目标回波
ML_value_refine_set_aver_weak = abs(ML_value_refine_set_aver - echo_strong_esti);

[value_row_max_refine, I_row_set_max_refine] = max(ML_value_refine_set_aver);
[~, I_col_max_refine] = max(value_row_max_refine);
I_row_max_refine = I_row_set_max_refine(I_col_max_refine);

figure
mesh(X,Y, 10*log10(abs(echo_strong_esti)./max(max(abs(echo_strong_esti)))));
ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
title('反演的强目标回波');

figure
mesh(X,Y, 10*log10(abs(ML_value_refine_set_aver_weak)./max(max(abs(echo_strong_esti)))));
ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
axis([-inf inf -inf inf -40 0]);
title('消除强目标回波后的残留结果');

%% =========== 8. 误差计算与切片对比 ===========
if exist('./result/para_base.mat', 'file')
    load('./result/para_base.mat', 'para_base');
    nor_delay_esti_base = para_base.nor_delay_esti_base;
    range_base_real = para_base.range_base_real;

    nor_delay_esti_wireless = (nor_delay_esti - nor_delay_esti_base);
    range_esti_wireless = nor_delay_esti_wireless/fs_up*c/2 + range_base_real/2;

    velocity_esti_wireless = nor_Doppler_esti*delta_f/24e9*c/2;
    nor_Doppler_real = 1.57;
    velocity_real_wireless = nor_Doppler_real*delta_f/24e9*c/2;
    del_velcity_esti = velocity_real_wireless - velocity_esti_wireless;

    range_real_wireless = 5.724;
    del_range_esti = (range_real_wireless/2 - range_esti_wireless);

    fprintf('误差结算：\n距离估计 = %.3f m (误差 %.3f)\n速度估计 = %.3f m/s (误差 %.3f)\n', ...
        range_esti_wireless, del_range_esti, velocity_esti_wireless, del_velcity_esti);
end

v_i_cal = [11.72];
delay_i_cal = [10.35];

figure
plot(Nor_delay_refineSearch_set, 10*log10(abs(ML_value_refine_set_aver(:,I_col_max_refine))./max(max(abs(echo_strong_esti)))), '-b*'); hold on;
plot(Nor_delay_refineSearch_set, 10*log10(abs(echo_strong_esti(:,I_col_max_refine))./max(max(abs(echo_strong_esti)))), '--ro');
xline(delay_i_cal, 'g--');
grid on; axis([-inf inf -30 2]);
set(gca,'FontSize',10);
xlabel('Normalized delay $\tau$', 'Interpreter', 'latex');
ylabel('Normalized average squared MF');
legend('Received','Reconstructed', 'True delay','Location','south');
annotation('textbox', [0.45 0.8 0.1 0.1], 'String', ['$\hat \tau = $', num2str(nor_delay_esti,'%10.5f')], ...
    'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 10,'Interpreter','latex');
title('时延切片对比');

figure
plot(Nor_Doppler_refineSearch_set, 10*log10(abs(ML_value_refine_set_aver(I_row_max_refine,:))./max(max(abs(echo_strong_esti)))), '-b*'); hold on;
plot(Nor_Doppler_refineSearch_set, 10*log10(abs(echo_strong_esti(I_row_max_refine,:))./max(max(abs(echo_strong_esti)))), '--ro');
xline(v_i_cal, 'g--');
grid on; axis([-inf inf -30 2]);
set(gca,'FontSize',10);
xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
ylabel('Normalized average squared MF');
legend('Received','Reconstructed', 'True doppler','Location','south');
annotation('textbox', [0.45 0.8 0.1 0.1], 'String', ['$\hat \nu = $', num2str(nor_Doppler_esti,'%10.5f')], ...
    'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 10,'Interpreter','latex');
title('多普勒切片对比');

disp('程序运行完毕。');
toc;