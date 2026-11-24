function [RxSignal, multiPathSignals, signalPower, noisePower] = f_time_vary_channel_v5(txSignal,...
    Nc, L, numPaths, channelCoefs, pathDelaysVector, pathNorDopplers, ...
    delayMax, sincFilterLen, SNR, addNoiseFlag)
% 多径衰落信道仿真函数
%
% 输入参数：
%   txSignal        - 发送信号（时域，上采样后）
%   Nc              - 子载波数量
%   L               - 过采样倍数
%   numPaths        - 多径数量(目标数)
%   channelCoefs    - 信道系数
%   pathDelaysVector- 各路径时延向量，给出每个采样点的时延，单位为上采样采样点
%   pathNorDopplers - 各路径归一化多普勒 nu，定义为 nu = f_d / delta_f
%   delayMax        - 最大时延
%   sincFilterLen   - 非整数时延滤波器长度
%   SNR             - 信噪比（dB）
%   addNoiseFlag    - 是否加噪
%
% 输出参数：
%   RxSignal         - 接收信号（含噪）
%   multiPathSignals - 各路径信号分量
%   signalPower      - 无噪声接收信号平均功率
%   noisePower       - 按 SNR 计算的复噪声平均功率
%
% Doppler 口径：
%   传入的 pathNorDopplers 是归一化 Doppler：nu = f_d / delta_f。
%   本函数在加 Doppler 相位时直接换算到上采样离散时间频率：
%       f_d / fs_up = (nu * delta_f) / (L * fs) = nu / (Nc * L)
%   因此相位为：
%       exp(j * 2*pi * (nu/(Nc*L)) * n)
%
% 修改历史:
%  ---------------------------------------------------------------------------------
%  版本      时间        修改者      修改内容                                  备注
%  ---------------------------------------------------------------------------------
%  v1    2026-01-26     [罗韬]     [初始版本,仅实现非时变时延和多普勒]
%  v2    2026-01-27     [罗韬]     [修改内容:
%                                  1.将sincFilterLen变量设置为可变参数
%                                  2.实现逐点添加时延，以适配时变时延]
%  v3    2026-02-03     [罗韬]     [修改内容:
%                                  1. 增加delayMax参数表示信道的最大时延
%                                  2. 实现过信道后的信号输入输出长度不一样
%  v4    2026-06-16     [罗韬]     [修改内容:
%                                  1. 增加是否加噪标志位
%  v5    2026-07-08     [Codex]   [修改内容:保持归一化 Doppler 输入，将 Doppler 单位转换整合到相位施加处]
%  v6    2026-11-24     [罗韬]     [修改内容:改用f_tv_delay_intfrac_v4，避免内存溢出]
%  ---------------------------------------------------------------------------------

    pathNorDopplers = pathNorDopplers(:);

    % 信号长度
    signalLen = length(txSignal);

    % 初始化输出
    RxSignal = zeros(signalLen + delayMax, 1);
    multiPathSignals = zeros(numPaths, signalLen + delayMax);

    % 上采样时域采样点索引。Doppler 相位中的频率单位是 cycles/sample。
    timeIdx = (0:signalLen - 1 + delayMax).';

    % 处理整数时延与分数时延
    intDelays = floor(pathDelaysVector(:, 1));
    fracDelays = pathDelaysVector - intDelays;
    addDe = kron(fracDelays(:, 1), ones(1, delayMax));
    maxDelays = [fracDelays, addDe];

    for p = 1:numPaths
        % 1. 先实现整数时延。
        IntDelayedSignal = circshift([txSignal; zeros(delayMax, 1)], intDelays(p, 1), 1);
        delaysCrishift(p, :) = circshift(maxDelays(p, :), intDelays(p, 1), 2);
        % delaysCrishift(p, :) = maxDelays(p, :);

        % 2. 逐点实现分数时延，以适配时变时延。
        % fracDelayedSignal = f_tv_delay_intfrac_v3(IntDelayedSignal, ...
        %     delaysCrishift(p, :), sincFilterLen);
        % v4版本占用更少内存
        fracDelayedSignal4 = f_tv_delay_intfrac_v4(IntDelayedSignal, ...
            delaysCrishift(p, :), sincFilterLen);

        % 3. 多普勒频移
        % pathNorDopplers(p) 是 nu = f_d / delta_f；
        % nu/(Nc*L) 是相对于上采样采样率 fs_up 的离散频率 f_d/fs_up。
        dopplerPhase = exp(1i * 2 * pi * (pathNorDopplers(p) / (Nc * L)) .* timeIdx);
        dopplerDelayedSignal = fracDelayedSignal .* dopplerPhase;

        % 4. 应用信道系数并累加。
        pathSignal = channelCoefs(p) * dopplerDelayedSignal;
        multiPathSignals(p, :) = pathSignal.';
        RxSignal = RxSignal + pathSignal;
    end

    % 5. 加噪
    signalPower = mean(abs(RxSignal).^2);
    SNRlinear = 10^(SNR / 10);
    noisePower = signalPower / SNRlinear;
    if addNoiseFlag
        noiseStd = sqrt(noisePower / 2);
        noise = noiseStd * (randn(size(RxSignal)) + 1i * randn(size(RxSignal)));
        RxSignal = RxSignal + noise;
    end

    % 转置 multiPathSignals 以保持输出维度一致
    multiPathSignals = multiPathSignals.';
end
