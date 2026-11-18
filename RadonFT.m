%% RadonFT.m
% 匀速目标 Radon-Fourier Transform（RFT）仿真与参数估计
%
% RFT 对每个候选速度执行两步处理：
%   1. Radon 部分：沿候选距离轨迹抽取各脉冲的距离压缩数据；
%   2. Fourier 部分：补偿候选轨迹对应的载频相位并相干求和。
% 候选速度等于真实速度时，目标包络和相位同时对齐，形成积累峰。
%
% 约定：径向速度 v>0 表示目标远离雷达，目标在 CPI 中心时刻的
% 距离记为 R0。搜索和回波仿真使用完全相同的慢时间参考。

close all;
clear;
% clc;
rng(2026);                         % [修改1] 固定随机种子，使结果可复现

%% 1. 单位及雷达参数
MHz = 1e6;
us  = 1e-6;
km  = 1e3;                         % [修改2] 1 km = 1e3 m，原代码误写为 10e3

c      = 3e8;                      % 光速，m/s
B      = 15*MHz;                   % LFM 带宽，Hz
fc     = 150*MHz;                  % 载频，Hz
lambda = c/fc;                     % [修改3] 显式定义波长，便于统一相位模型
fs     = 40*MHz;                   % 快时间采样率，Hz
Tp     = 50*us;                    % 脉冲宽度，s
Ts     = 1/fs;                     % 快时间采样间隔，s
Kr     = B/Tp;                     % LFM 调频率，Hz/s

pulseNum = 512;                    % 脉冲数
Tcpi     = 1;                      % 相干处理时间，s
PRI      = Tcpi/pulseNum;          % 脉冲重复间隔，s
PRF      = 1/PRI;                  % 脉冲重复频率，Hz

pulseSampleNum = round(Tp*fs);
recordSampleNum = 4096;            % [修改4] 按 2 km 场景设置足够的快时间记录长度
guardSampleNum  = 128;
fastTime = (-pulseSampleNum/2:pulseSampleNum/2-1)*Ts;

% [修改5] CPI 中心作为 t=0；仿真和 RFT 搜索均使用 slowTime，
% 消除原代码 (k-PulseNum/2)*PRI 与 k*PRI 混用造成的距离偏差。
slowTime = ((0:pulseNum-1)-(pulseNum-1)/2)*PRI;
rangeBinSpacing = c*Ts/2;

%% 2. 三个匀速目标
R0True = 2*km;                      % CPI 中心时刻距离，m
vTrue  = [300];           % 径向速度，m/s；正值表示远离
targetNum = numel(vTrue);

rangeTrue = R0True + vTrue(:)*slowTime;
delaySample = round(2*rangeTrue/(c*Ts));

% [修改6] 在生成数组前检查所有目标回波是否位于记录范围内。
firstEchoSample = guardSampleNum + min(delaySample(:)) + 1;
lastEchoSample  = guardSampleNum + max(delaySample(:)) + pulseSampleNum;
if firstEchoSample < 1 || lastEchoSample > recordSampleNum
    error('目标回波超出快时间记录，请增大 recordSampleNum。');
end

%% 3. 生成复基带 LFM 回波
% [修改7] 删除未被后续处理使用、且会产生带通欠采样歧义的射频波形；
% 直接生成复基带 LFM 信号，其实部和虚部分别对应 I、Q 两路。
txChirp = exp(1j*pi*Kr*fastTime.^2);
echo = complex(zeros(pulseNum,recordSampleNum));

for p = 1:pulseNum
    for q = 1:targetNum
        firstSample = guardSampleNum + delaySample(q,p) + 1;
        echoIndex = firstSample:firstSample+pulseSampleNum-1;

        % [修改8] 使用标准窄带“停走停”模型。回波相位采用负号，
        % 后面的 RFT 使用正号补偿；二者必须成对一致。
        targetPhase = exp(-1j*4*pi*rangeTrue(q,p)/lambda);
        echo(p,echoIndex) = echo(p,echoIndex) + txChirp*targetPhase;
    end
end

SNRdB = 10;
% [修改9] SNR 以单个单位幅度目标的脉内功率为基准，不再用包含大量零采样的整个回波矩阵测量功率，避免噪声功率被低估。
noisePower = 10^(-SNRdB/10);
noise = sqrt(noisePower/2).* (randn(size(echo)) + 1j*randn(size(echo)));
echoNoisy = echo + noise;

%% 4. 线性脉冲压缩及距离标定
% [修改10] 使用标准时间反转共轭匹配滤波器和线性卷积，避免原代码
% 多次 FFT 写法中的循环卷积回绕及峰值索引歧义。
matchedFilter = conj(fliplr(txChirp));
pulseCompressed = complex(zeros(pulseNum,recordSampleNum));
for p = 1:pulseNum
    pulseCompressed(p,:) = conv(echoNoisy(p,:),matchedFilter,'same');
end

% [修改11] 用零时延标定脉冲确定 conv(...,'same') 的峰值索引，再建立具有物理单位的距离轴。
calibrationEcho = complex(zeros(1,recordSampleNum));
calibrationIndex = guardSampleNum+1:guardSampleNum+pulseSampleNum;
calibrationEcho(calibrationIndex) = txChirp;
calibrationCompressed = conv(calibrationEcho,matchedFilter,'same');
[~,zeroDelayPeakBin] = max(abs(calibrationCompressed));
rangeAxis = ((1:recordSampleNum)-zeroDelayPeakBin)*rangeBinSpacing;

figure('Name','脉冲压缩结果');
plot(rangeAxis/km,abs(pulseCompressed(1,:)),'k');
hold on;
plot(rangeAxis/km,abs(pulseCompressed(end,:)),'b');
xlim([1.6 2.4]);
xlabel('距离（km）');
ylabel('幅值');
title('第 1 个与最后一个脉冲的脉冲压缩结果');
legend('第 1 个脉冲','最后一个脉冲');
grid on;

figure('Name','距离-慢时间数据');
imagesc(rangeAxis/km,slowTime,abs(pulseCompressed));
axis xy;
xlim([1.6 2.4]);
xlabel('距离（km）');
ylabel('慢时间（s）');
title('脉冲压缩后的距离走动');
colorbar;

%% 5. 常规 MTD（用于与 RFT 对比）
% [修改12] 明确沿第 1 维（慢时间维）进行 FFT 和 fftshift。
mtd = fftshift(fft(pulseCompressed,[],1),1);
dopplerAxis = (-pulseNum/2:pulseNum/2-1)*(PRF/pulseNum);
displayRangeBins = find(rangeAxis >= 1.6*km & rangeAxis <= 2.4*km);

figure('Name','常规 MTD 结果');
imagesc(rangeAxis(displayRangeBins)/km,dopplerAxis, ...
    20*log10(abs(mtd(:,displayRangeBins))/max(abs(mtd(:)))+eps));
axis xy;
xlabel('距离（km）');
ylabel('多普勒频率（Hz）');
title('未校正距离走动时的 MTD 结果（dB）');
colorbar;

%% 6. 匀速目标 RFT 参数搜索
deltaV = 1;                         % 速度搜索步长，m/s
velocitySearch = -500:deltaV:500;
rangeSearchMin = 1.6*km;
rangeSearchMax = 2.4*km;

% [修改16] 将速度搜索、Radon 轨迹抽取、Fourier 相位补偿和相干积累
% 封装到独立的 RFT.m 函数中。主脚本只负责准备输入和处理输出。
[rftResult,rftRangeAxis] = RFT(pulseCompressed,rangeAxis,slowTime, ...
    velocitySearch,lambda,[rangeSearchMin rangeSearchMax]);

rftPower = abs(rftResult).^2;
rftPowerNormalized = rftPower/max(rftPower(:));

%% 7. 三目标参数估计
% [修改15] 原代码只显示 RFT 曲面，没有读取参数。这里采用二维非极大值
% 抑制依次提取三个峰，并避免把同一目标的邻近旁瓣重复识别为多个目标。
peakCount = targetNum;
workingPower = rftPower;
vEstimate = zeros(1,peakCount);
R0Estimate = zeros(1,peakCount);
peakPower = zeros(1,peakCount);
velocityGuardBins = max(5,round(5/deltaV));
rangeGuardBins = max(4,ceil((c/(2*B))/rangeBinSpacing));

for q = 1:peakCount
    [peakPower(q),peakLinearIndex] = max(workingPower(:));
    [velocityIndex,localRangeIndex] = ind2sub(size(workingPower),peakLinearIndex);
    vEstimate(q) = velocitySearch(velocityIndex);
    R0Estimate(q) = rftRangeAxis(localRangeIndex);

    velocitySuppress = max(1,velocityIndex-velocityGuardBins) ...
                     : min(size(workingPower,1),velocityIndex+velocityGuardBins);
    rangeSuppress = max(1,localRangeIndex-rangeGuardBins) ...
                  : min(size(workingPower,2),localRangeIndex+rangeGuardBins);
    workingPower(velocitySuppress,rangeSuppress) = 0;
end

% 按估计速度排序，便于与真实目标一一比较。
[vEstimate,sortIndex] = sort(vEstimate);
R0Estimate = R0Estimate(sortIndex);
peakPower = peakPower(sortIndex);

fprintf('\n========== RFT 三目标参数估计结果 ==========\n');
fprintf('目标    真实距离/km   估计距离/km   真实速度/(m/s)   估计速度/(m/s)\n');
for q = 1:peakCount
    fprintf('%2d        %7.4f        %7.4f          %8.2f          %8.2f\n', ...
        q,R0True/km,R0Estimate(q)/km,vTrue(q),vEstimate(q));
end

figure('Name','RFT 参数搜索平面');
imagesc(rftRangeAxis/km,velocitySearch, 10*log10(rftPowerNormalized+eps));
axis xy;
xlabel('CPI 中心时刻距离（km）');
ylabel('搜索速度（m/s）');
title('RFT 归一化相干积累功率（dB）');
colorbar;
hold on;
plot(R0Estimate/km,vEstimate,'rx','MarkerSize',11,'LineWidth',2);

figure('Name','RFT 三维积累结果');
mesh(rftRangeAxis/km,velocitySearch,rftPowerNormalized);
xlabel('CPI 中心时刻距离（km）');
ylabel('搜索速度（m/s）');
zlabel('归一化积累功率');
title('匀速目标 RFT 积累结果');
grid on;
view(45,35);

% 提示：载频相位对应的单独无模糊速度范围为 ±lambda*PRF/4。
% RFT 同时利用距离走动信息，能够在一定条件下进一步区分多普勒模糊速度。
