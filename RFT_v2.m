function [rftResult, rftRangeAxis] = RFT_v2(pulseCompressed, rangeAxis, slowTime, velocitySearch, lambda, rangeSearchLimits)
% RFT_V2 逐速度执行 RFT，用距离维 FFT/IFFT 实现分数格轨迹取样。
% pulseCompressed：复数脉压数据，行为脉冲，列为距离采样点。
% rangeAxis：等间隔递增的距离轴（m，至少两点）。
% slowTime：各脉冲相对 t=0 的时刻（s）。
% velocitySearch：候选径向速度（m/s，远离为正）。
% lambda：波长（m）。
% rangeSearchLimits：[最小距离 最大距离]（m），省略或传 [] 时搜索全轴。
% rftResult：行对应候选速度，列对应 rftRangeAxis 上 t=0 的候选距离。

    rangeAxis = rangeAxis(:).';
    slowTime = slowTime(:);
    rangeStep = median(diff(rangeAxis));
    rangeBinCount = size(pulseCompressed, 2);

    % 搜索窗口只限制 t=0 的候选距离。轨迹取样仍可访问整个距离记录。
    rangeBins = 1:rangeBinCount;
    if nargin >= 6 && ~isempty(rangeSearchLimits)
        rangeBins = find(rangeAxis >= rangeSearchLimits(1) & rangeAxis <= rangeSearchLimits(2));
    end
    rftRangeAxis = rangeAxis(rangeBins);
    rftResult = complex(zeros(numel(velocitySearch), numel(rangeBins)));
    if isempty(rftResult)
        return;
    end

    % 每行沿距离维做一次 FFT，结果由所有候选速度共用。
    % FFT 时移是周期性的。在记录末尾补零，为正、负时移留出余量。
    maxShiftBins = ceil(max(abs(velocitySearch)) * max(abs(slowTime)) / rangeStep);
    fftLength = rangeBinCount + 2 * maxShiftBins;
    rangeSpectrum = fft(pulseCompressed, fftLength, 2);
    % 频率单位为 cycles/sample，正负频率按 MATLAB fft 的输出顺序排列。
    frequency = ifftshift(-floor(fftLength/2):ceil(fftLength/2)-1) / fftLength;

    for k = 1:numel(velocitySearch)
        rangeChange = velocitySearch(k) * slowTime;
        binShift = rangeChange / rangeStep;  % 位移以采样格为单位，保留小数部分。

        % Radon：第 m 行在 r0 + v*t_m 取值。若位移为 d 格，
        % 频域乘 exp(+j*2*pi*f*d) 后做 IFFT，得到 x[n+d] 在整数 n 上的值。
        % 例如 d=0.4 时，输出第 n 格就对应原距离剖面的 n+0.4 格。
        shiftPhase = exp(2j*pi * (binShift * frequency));
        shiftedProfiles = ifft(rangeSpectrum .* shiftPhase, [], 2);
        trajectoryValues = shiftedProfiles(:, rangeBins);

        % 轨迹超出实际记录的样本按零处理，避免使用周期延拓出来的值。
        queryBins = rangeBins + binShift;
        trajectoryValues(queryBins < 1 | queryBins > rangeBinCount) = 0;

        % Fourier：按原 RFT 约定，用未取整的位移补偿载波相位，再沿脉冲求和。
        % 这里使用非共轭转置，保留正号补偿 exp(+j*4*pi*v*t/lambda)。
        phaseCompensation = exp(4j*pi * rangeChange / lambda);
        rftResult(k, :) = phaseCompensation.' * trajectoryValues;
    end
end
