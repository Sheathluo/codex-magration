%%  版本信息
%% 2025/12/23 V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2025  BUAA
%======================================================
% Designer：Ni Yuanhan
% Design Date：20251223
% Modifier：Ni Yuanhan
% Modifying Date：20251223
% Function Description：模拟感知信道，非整数时延、非整数多普勒
%======================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
% rng(1); % 固定随机数种子，确保每次加噪序列完全一致，方便算法排查
tic

% 检查并创建结果保存目录
if ~exist('./result', 'dir')
    mkdir('./result');
end

%% *********** 1. 读取参数 *************************
para_path = './data/para_set_pulsedAFDM_N128_ISAC_time_20251124_X410.mat';
if ~exist(para_path, 'file')
    error('未找到参数文件，请确认路径: %s', para_path);
end
load(para_path, 'para_set');

L           = para_set.L;
beta        = para_set.beta;
L_delay     = para_set.L_delay;
L_Doppler   = para_set.L_Doppler;
flag_pulse  = para_set.flag_pulse;
flag_RRC    = para_set.flag_RRC;
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
N_preSuffix = para_set.N_preSuffix; % 波形添加的前缀或后缀长度
N_tapFiler  = para_set.N_tapFiler;  % 收发等效脉冲成形滤波器半边抽头数量
delta_f     = para_set.delta_f;

%% *********** 2. 设置感知信道参数 *************************
Num_target = 1;
chan_coef = [1].';
v_i_cal = [11.72];
delay_i_cal = [10.35];   % 固定时延（归一化）

% 距离徙动，如果v=0，则不考虑徙动，如果不为零则考虑徙动
velocity = [0, 0];          % 速度
% velocity = v_i_cal .* delta_f * c / fc;

% 速度引起的时延变化 (双程)，徙动在采样点之间变化
delta_tau = (1:Nsym) * 2 .* (velocity' .* L * (Nc + Ncp + 2 * N_preSuffix) / c);
% delta_tau = (1:(Nsym * L * (Nc + Ncp + 2 * N_preSuffix))) * 2 .* (velocity' / c);

% 考虑时延在采样点间变化
tau = delay_i_cal' + delta_tau;
% 考虑符号内时延不变，符号间时延变化
tau_vector_sym = kron(tau, ones(1, L * (Nc + Ncp + 2 * N_preSuffix)));

% 处理两头的数据
tau_vector_head = tau(:,1) .* ones(1, N_preSuffix);
tau_vector_end = tau(:,end) .* ones(1, N_preSuffix);
tau_vector = [tau_vector_head, tau_vector_sym, tau_vector_end];
delay_max = Ncp * L - 1;

SNR_in_dB = 10; % 输入波形信噪比

% range_migration_para = struct();
% range_migration_para.tau = tau;
% range_migration_para.delta_tau_migration = delta_tau;
% save('./data/para_range_migration_vactor.mat', 'range_migration_para');
%% *********** 3. 读取发送信号 *************************
filenamere = './data/Tx_pulsedAFDM_N128_ISAC_time_20251124_X410_re.dat';
filenameim = './data/Tx_pulsedAFDM_N128_ISAC_time_20251124_X410_im.dat';

fidre = fopen(filenamere, 'r');
if fidre > 0
    Rx_Sigre = fread(fidre, 'float32');
    fclose(fidre);
else
    error('无法读取发送信号实部文件: %s', filenamere);
end

fidim = fopen(filenameim, 'r');
if fidim > 0
    Rx_Sigim = fread(fidim, 'float32');
    fclose(fidim);
else
    error('无法读取发送信号虚部文件: %s', filenameim);
end

% 对齐实部与虚部长度
len_data = min(length(Rx_Sigre), length(Rx_Sigim));
Tx_signal = complex(Rx_Sigre(1:len_data), Rx_Sigim(1:len_data));

%% *********** 4. 经过多径信道 ******************
% 添加时变时延
pathDelaysVector = tau_vector;
sincFilterLen = 500;
[Rx_signal_AFDM, ~] = f_time_vary_channel_v3(Tx_signal, Nc, L, Num_target, ...
    chan_coef, pathDelaysVector, v_i_cal, delay_max, sincFilterLen, 1);

%% *********** 5. 加噪 ***********************
P_rx_cal = mean(abs(Rx_signal_AFDM).^2);
P_noise = P_rx_cal / (10^(SNR_in_dB/10));
Len = length(Rx_signal_AFDM);

if flag_RRC == 1
    noise_AFDM = sqrt(L) * sqrt(P_noise/2) * (randn(Len,1) + 1i*randn(Len,1));
else
    noise_AFDM = sqrt(P_noise/2) * (randn(Len,1) + 1i*randn(Len,1));
end

% 加噪
% Rx_signal_AFDM = Rx_signal_AFDM + noise_AFDM;

%% *********** 6. 保存接收信号 ***************************
Rx_re = real(Rx_signal_AFDM);
Rx_im = imag(Rx_signal_AFDM);

out_filenamere = './result/Rx_simSensing_pulsedAFDM_N128_ISAC_time_20251124_X410_re.dat';
out_filenameim = './result/Rx_simSensing_pulsedAFDM_N128_ISAC_time_20251124_X410_im.dat';

fidre_out = fopen(out_filenamere, 'w');
if fidre_out > 0
    fwrite(fidre_out, Rx_re, 'float32');
    fclose(fidre_out);
else
    error('接收信号实部写入失败，请检查 ./result/ 路径权限！');
end

fidim_out = fopen(out_filenameim, 'w');
if fidim_out > 0
    fwrite(fidim_out, Rx_im, 'float32');
    fclose(fidim_out);
else
    error('接收信号虚部写入失败，请检查 ./result/ 路径权限！');
end

disp('感知信道模拟完成，接收信号已成功导出至 result 文件夹。');
toc;