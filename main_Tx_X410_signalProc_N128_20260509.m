%%  版本信息
%%
%% 2026/05/09 V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026  BUAA
%======================================================
% Designer：Ni Yuanhan
% Design Date：20250626
% Modifier：Ni Yuanhan
% Modifying Date：20250805
% Function Description：对比OFDM和AFDM在强、弱目标条件下的测速误差。
% Modification：v11:在v10的基础上，优化代码，适配X410等软件无线电实验
%               v12:在v11的基础上，增大采样率，适配X410等软件无线电实验
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
rng(1); % 【新增】：固定随机数种子，确保每次生成的发送数据和波形百分之百一致
tic

% 检查并创建保存数据的目录，防止由于缺少文件夹导致报错
if ~exist('./data', 'dir')
    mkdir('./data');
end

%% 脉冲成形参数
L = 2;                          % 过采样倍数
beta = 0.35;                    % 滚降因子
L_delay = 2;                    % 时延过采样倍数
L_Doppler = 2;                  % 多普勒轴过采样倍数
flag_pulse = 1;
flag_RRC = 1;
N_preSuffix = 20;               % 波形添加的前缀或后缀长度
N_tapFiler = 20;                % 收发等效脉冲成形滤波器半边抽头数量 2*N_tapFiler*L + 1
% fs_up = 49.152e6/8;
fs_up = 80e6/8;
fs = fs_up/L;

%% 通信（AFDM）参数
c = 3e8;                        % 光速
fc = 24*1e9;                    % Hz,载波频率
lamda = c/fc;                   % 载波波长

Nc_0 = 128;
Nc = Nc_0*1;                    % 子载波个数
Nsym = 500;                     % 码元个数/帧数目
Mod = 16;                       % 16-QAM调制
B = fs;
delta_f = B/Nc;
T = 1/delta_f;                  % s, 码元持续时间
Ncp = 20;                       % s, 循环前缀点数
Tcp = Ncp/B;                    % s, 循环前缀时间
T_AFDM = T+Tcp;                 % s, 整个AFDM符号持续时间

N_delay = Nc;
N_Doppler = Nc;
N = Nc;
M = Nc;

%% 雷达 & 目标参数
Num_target = 2;                 % 目标（多径）数目
delay_max = Num_target;
alpha_max = 2;                  % 被delta_f归一化的多普勒频率
xi_v = 4;                       % 由于非整数倍归一化多普勒，造成等效信道宽度扩展
Q = (delay_max + 1)*(2*(alpha_max + xi_v)+1) - 1;

c1 = (8)/(2*Nc);
c2 = 0;

doppler_max = alpha_max*delta_f;
doppler_max_perPath = [Nc/(Nc+Ncp)*delta_f, Nc/(Nc+Ncp)*delta_f, doppler_max];
delta_tau = T/Nc; % 时延分辨率
delta_fd = 1/(T_AFDM*Nsym); % 多普勒分辨率
delta_range_com = delta_tau*c; % 距离分辨率（通信，单程）
delta_velocity_com = delta_fd/fc*c; % 速度分辨率（通信，单程）m/s
delta_range_rad = delta_tau*c/2; % 距离分辨率（雷达，双程）
delta_velocity_rad = delta_fd/fc*c/2; % 速度分辨率（雷达，双程）m/s

Max_unamb_range = c*T/2; % 最大不模糊距离（雷达，双程）
Max_ISI_free_range = c*Tcp/2; % 最大无ISI距离（雷达，双程）,时延不超过Tcp
Max_unamb_velocity = c/(2*fc*T_AFDM); % 最大不模糊速度
Max_ICI_free_velocity = B*c/(20*Nc*fc);

velocity_delta_f = delta_f/(2*fc)*c;
velocity_1_5_delta_f = 1.5*delta_f/(2*fc)*c;
velocity_2_delta_f = 2*delta_f/(2*fc)*c;

% SNR与强弱目标参数
SNR_input_dB = 25; % 接收到信号信噪比 dB
SNR_input_set = 10.^(SNR_input_dB/10);
SNR_Pilot_dB = 40; % 导频信噪比
PSWR_dB = 21; % 强弱目标功率比
PSWR_set = 10.^(PSWR_dB/10);
SNR_weak_dB = 20;
SNR_weak_set = 10.^(SNR_weak_dB/10);
P_noise = 1;

%% 生成成形滤波器响应参数
mu_4 = 1.32; % 16-QAM
delay_nor_set = 0:1/L_delay:N*L-1/L_delay;
Doppler_nor_set = -N*L/2:1/L_Doppler:N*L/2-1/L_Doppler;

%% 通信/雷达发射端 - 产生基带信号
bit_num = Nc;
JSR_dB = 10;
Datain = randi([0 1], bit_num * Nsym * log2(Mod),1);
ModPara = struct('Mod',Mod,'Tb',1/20/1e6,'fs',20*1e6,'Rb',20*1e6);

Datain_res = reshape(Datain,log2(Mod),[]).';
data_temp = bi2de(Datain_res);  % 二进制转十进制
BPSKsig = (qammod(data_temp,Mod,'UnitAveragePower',true)).';  % QAM调制
kron_vec = ones(1,round(ModPara.fs * ModPara.Tb));
BPSKsig_upsamp = kron(BPSKsig,kron_vec);

BPSKsig_data_reshape = reshape(BPSKsig_upsamp,length(BPSKsig_upsamp)/Nsym,[]);
BPSKsig_reshape = BPSKsig_data_reshape;

% 生成单载波 SC
Tx_signal_SC = reshape(BPSKsig_reshape,[],1);

% 生成全CP AFDM
waveform_index = 'AFDM';
N_valid = Nc;
vec_BPSKsig_reshape = BPSKsig_reshape;
N_AFDM = N;
M_AFDM = N_AFDM;
Nsym_AFDM = Nsym;
[~, transmit_data_noCP_AFDM, transmit_data_addCP_AFDM, ~] = f_waveform_generator_ADDM_v3(waveform_index, c1, c2, M_AFDM, N_AFDM, Ncp,Nsym_AFDM, vec_BPSKsig_reshape);
Tx_signal_AFDM = transmit_data_noCP_AFDM;

% transmit_data_noCP_AFDM_ser = reshape(transmit_data_noCP_AFDM, [], 1);
% Nfft = length(transmit_data_noCP_AFDM_ser);
% specAFDM = fftshift(fft(transmit_data_noCP_AFDM_ser, Nfft)).';
% f = (-Nfft/2:Nfft/2-1)*(fs/Nfft);
% figure;
% plot(f,10*log10(abs(specAFDM)));

%% *********** 脉冲成形 方式三：先加CP，再加前缀和后缀，再成形 *****************************
transmit_data_CP_preSuffix_AFDM = [transmit_data_noCP_AFDM(end-Ncp-N_preSuffix+1:end-Ncp,:); transmit_data_noCP_AFDM(end-Ncp+1:end,:); transmit_data_noCP_AFDM; transmit_data_noCP_AFDM(1:N_preSuffix,:)];

% 符号过采样
x_sigCP_preSuffix = reshape(transmit_data_CP_preSuffix_AFDM,[],1).';
x_upSamCP_preSuffix = reshape([sqrt(L)*x_sigCP_preSuffix; zeros(L-1,(Nc+Ncp+2*N_preSuffix)*Nsym)], [],1);

X_upSamCP_preSuffix = 1/sqrt(length(x_upSamCP_preSuffix))*abs(fft(x_upSamCP_preSuffix));

norm_ps_TxRx = zeros(1,1);
g_ps_TxRx = zeros(Nc*L+1,1);
g_ps_TxRx_shift = zeros(Nc*L,1);
N_tapFiler_tx = 0;

% 生成成形滤波器响应
if flag_pulse == 1
    if flag_RRC == 1
        % 直接用于线卷积滤波的 脉冲成形响应
        N_tapFiler_tx = N_tapFiler/2;
        g_ps = rcosdesign(beta,N_tapFiler,L,'sqrt').'; % 根升余弦
        g_ps = g_ps/norm(g_ps);

        % 等效循环卷积的 收发等效 脉冲成形响应
        conv_filter = conv(conj(flip(g_ps)),g_ps);
        norm_ps_TxRx = norm(conv_filter);

        % 越界保护与自适应长度填充
        start_idx = Nc*L/2 - N_preSuffix*L + 1;
        end_idx = start_idx + length(conv_filter) - 1;
        if end_idx <= length(g_ps_TxRx)
            g_ps_TxRx(start_idx:end_idx) = conv_filter;
        else
            error('滤波器抽头超出分配的信号长度，请适当增大 Nc 或减小 N_preSuffix/N_tapFiler。');
        end
        g_ps_TxRx_shift = circshift(g_ps_TxRx(1:Nc*L,1),-Nc*L/2);
    else
        N_tapFiler_tx = N_tapFiler;
        g_ps = rcosdesign(beta,2*N_tapFiler,L,'normal').'; % 升余弦
        g_ps = g_ps/norm(g_ps);
        norm_ps_TxRx = norm(g_ps);

        start_idx = Nc*L/2 - N_preSuffix*L + 1;
        end_idx = start_idx + length(g_ps) - 1;
        if end_idx <= length(g_ps_TxRx)
            g_ps_TxRx(start_idx:end_idx) = g_ps;
        else
            error('滤波器抽头超出分配的信号长度，请检查参数。');
        end
        g_ps_TxRx_shift = circshift(g_ps_TxRx(1:Nc*L,1),-Nc*L/2);
    end
else
    g_ps = [zeros(N_tapFiler*L,1); 1; zeros(N_tapFiler*L,1)];
    g_ps_TxRx = g_ps;
end

P_ps_TxRx = zeros(Nc*L,Nc*L);
for i = 0:Nc*L-1
    P_ps_TxRx(:,i+1) = circshift(g_ps_TxRx_shift,i);
end

% 利用线性卷积，实现成形滤波
x_pulseShaped_CP_preSuffix = conv(g_ps,x_upSamCP_preSuffix);
% 成形滤波后的波形频谱
X_pulseShaped_CP_preSuffix = 1/sqrt(length(x_pulseShaped_CP_preSuffix))*abs(fft(x_pulseShaped_CP_preSuffix));

%% ************ 等效的脉冲成形信号（用于匹配滤波的参考信号） ************
p_rrc_trunc = zeros(Nc*L,1);
start_idx_rrc = Nc*L/2 - N_preSuffix*L/2 + 1;
end_idx_rrc = start_idx_rrc + length(g_ps) - 1;
p_rrc_trunc(start_idx_rrc:end_idx_rrc) = g_ps;

% 循环移位
p_rrc = circshift(p_rrc_trunc, -Nc*L/2);
P_rrc = zeros(Nc*L,Nc*L);
for i = 1:Nc*L
    P_rrc(:,i) = circshift(p_rrc,i-1);
end

x_upSam = kron(sqrt(L)*transmit_data_noCP_AFDM, [1;zeros(L-1,1)]);
% 经过脉冲成形后的发送波形
Tx_signal_pulsedAFDM_eff = P_rrc * x_upSam;

%% ************* 保存参数、画图与生成 dat 文件 ****************
para_set = struct();
para_set.Nc = Nc;
para_set.L = L;
para_set.beta = beta;
para_set.L_delay = L_delay;
para_set.L_Doppler = L_Doppler;
para_set.flag_pulse = flag_pulse;
para_set.flag_RRC = flag_RRC;
para_set.N_preSuffix = N_preSuffix;
para_set.N_tapFiler = N_tapFiler;
para_set.fs_up = fs_up;
para_set.fs = fs;
para_set.c = c;
para_set.fc = fc;
para_set.Nsym = Nsym;
para_set.Mod = Mod;
para_set.B = B;
para_set.Ncp = Ncp;
para_set.c1 = c1;
para_set.c2 = c2;
para_set.delta_f = delta_f;
save('./data/para_set_pulsedAFDM_N128_ISAC_time_20251124_X410.mat', 'para_set');

% 保存发射待信号 (归一化)
x_pulseShaped_CP_preSuffix_nor = 0.1 * x_pulseShaped_CP_preSuffix ./ max(abs(x_pulseShaped_CP_preSuffix));

figure;
plot(20*log10(abs(x_pulseShaped_CP_preSuffix_nor)));
title('发送信号-频谱');

figure;
plot(abs(x_pulseShaped_CP_preSuffix_nor));
title('发送信号-时域-包络');

% 分离实部与虚部
Tx_re = real(x_pulseShaped_CP_preSuffix_nor);
Tx_im = imag(x_pulseShaped_CP_preSuffix_nor);

% ======== 安全写入发送信号 ========
filenamere = './data/Tx_pulsedAFDM_N128_ISAC_time_20251124_X410_re.dat';
fidre = fopen(filenamere,'w');
if fidre > 0
    fwrite(fidre, Tx_re, 'float32');
    fclose(fidre);
else
    error('发送信号实部地址写入失败，请检查路径权限！');
end

filenameim = './data/Tx_pulsedAFDM_N128_ISAC_time_20251124_X410_im.dat';
fidim = fopen(filenameim,'w');
if fidim > 0
    fwrite(fidim, Tx_im, 'float32');
    fclose(fidim);
else
    error('发送信号虚部地址写入失败，请检查路径权限！');
end

% ======== 安全写入参考信号 ========
Tx_signal_pulsedAFDM_eff_vec = reshape(Tx_signal_pulsedAFDM_eff, [], 1);
Tx_re_ref = real(Tx_signal_pulsedAFDM_eff_vec);
Tx_im_ref = imag(Tx_signal_pulsedAFDM_eff_vec);

filenamere_ref = './data/Tx_ConfencedSignal_pulsedAFDM_N128_ISAC_time_20251124_X410_re.dat';
fidre_ref = fopen(filenamere_ref,'w');
if fidre_ref > 0
    fwrite(fidre_ref, Tx_re_ref, 'float32');
    fclose(fidre_ref);
else
    error('参考信号实部地址写入失败，请检查路径权限！');
end

filenameim_ref = './data/Tx_ConfencedSignal_pulsedAFDM_N128_ISAC_time_20251124_X410_im.dat';
fidim_ref = fopen(filenameim_ref,'w');
if fidim_ref > 0
    fwrite(fidim_ref, Tx_im_ref, 'float32');
    fclose(fidim_ref);
else
    error('参考信号虚部地址写入失败，请检查路径权限！');
end

disp('脚本执行完毕，波形已成功导出。');
toc;