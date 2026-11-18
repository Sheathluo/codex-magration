function [rftResult, rftRangeAxis] = RFT(pulseCompressed, rangeAxis, slowTime, velocitySearch, lambda, rangeSearchLimits)
% RFT 对匀速目标执行标准离散 Radon-Fourier 变换（最近距离单元取样）。
% pulseCompressed：复数距离压缩数据，行为脉冲，列为距离单元。
% rangeAxis：等间隔递增的距离轴（m，至少两点）；
% slowTime：各脉冲的时刻（s）。
% velocitySearch：候选径向速度（m/s，远离为正）；
% lambda：波长（m）。
% rangeSearchLimits：[最小距离 最大距离]（m），省略或传 [] 时搜索整条距离轴。
% rftResult 的行为候选速度，列为 rftRangeAxis 上 t=0 时的候选距离。

    rangeAxis = rangeAxis(:).';
    slowTime = slowTime(:);
    [pulseCount, rangeBinCount] = size(pulseCompressed);
    rangeStep = median(diff(rangeAxis));  % 距离单元间隔（m）

    % 搜索范围限制的是 t=0 的距离；轨迹取样仍可使用整条距离轴。
    rangeBins = 1:rangeBinCount;
    if nargin >= 6 && ~isempty(rangeSearchLimits)
        rangeBins = find(rangeAxis >= rangeSearchLimits(1) & rangeAxis <= rangeSearchLimits(2));
    end
    rftRangeAxis = rangeAxis(rangeBins);
    rftResult = complex(zeros(numel(velocitySearch), numel(rangeBins)));

    % 每个候选速度对应输出的一行；这一行同时搜索所有候选初始距离。
    for k = 1:numel(velocitySearch)
        rangeChange = velocitySearch(k) * slowTime;  % 各脉冲相对 t=0 的距离变化 v*t
        binShift = round(rangeChange / rangeStep);   % 换算成整数距离单元

        % 回波相位约定为 exp(-j*4*pi*R/lambda)，用未取整的距离变化补偿相位。
        phaseCompensation = exp(1j * 4*pi * rangeChange / lambda);

        for m = 1:pulseCount
            % Radon：从第 m 个脉冲沿 R0+v*t 取样。同一速度下，各 R0 平移相同格数。
            sampleBins = rangeBins + binShift(m);
            valid = sampleBins >= 1 & sampleBins <= rangeBinCount;

            % Fourier：补偿后直接累加复数回波。越界位置跳过，相当于补零。
            echo = pulseCompressed(m, sampleBins(valid));
            rftResult(k, valid) = rftResult(k, valid) + echo * phaseCompensation(m);
        end
    end
end
