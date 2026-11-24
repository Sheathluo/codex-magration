function y = f_tv_delay_intfrac_v4(x, tau, M)
%------------------------------------------------------------
%------------------------------------------------------------
% 时变时延（整数 + 分数），每个采样点对应不同的分数延时滤波器windowed-sinc 实现
%
% x   : 输入信号 [N x 1] (实或复)
% tau : 总时延轨迹 [N x 1]，任意实数（可>1，可<0）
% M   : sinc 半长度，核长 = 2M+1
%
% y   : 输出信号 [N x 1]
%
% 模型：
%   tau[n] = D[n] + delta[n], delta[n]∈[0,1)
%   y[n] = sum_{m=-M..M} x[n - D[n] - m] * sinc(m-delta[n]) * w[m]
%------------------------------------------------------------
% 修改历史:
%  ---------------------------------------------------------------------------------
%  版本      时间        修改者      修改内容                                  备注
%  ---------------------------------------------------------------------------------
%  v3    2026-05-20     [罗韬]     [极致矢量化优化的时变分数时延滤波器]
%  v4    2026-11-24     [罗韬]     [将全帧矩阵改为分块矩阵，避免长信号内存不足。]
%
%------------------------------------------------------------

    x = x(:);
    tau = tau(:);
    N = length(x);
    y = zeros(N, 1, 'like', x);

    % ===== Step 1: 按核长确定分块大小 =====
    % 每个二维临时矩阵最多约 2^18 个元素；double 实矩阵约 2 MiB，
    % 复矩阵约 4 MiB。该值限制单个矩阵，而非所有临时变量的总内存。
    % 若核长本身超过此值，则每块只计算一个输出采样点。
    maxMatrixElements = 2^18;
    kernelLength = 2*M + 1;
    blockSize = max(1, floor(maxMatrixElements / kernelLength));
    m = (-M:M).';

    % ===== Step 2: 逐块生成滤波核并计算输出 =====
    for first = 1:blockSize:N
        last = min(first + blockSize - 1, N);
        n = first:last;             % 全帧输出索引，不能从 1 重新计数

        if isscalar(tau)
            tau_block = tau;
        else
            tau_block = tau(first:last).';
        end
        D = round(tau_block);
        delta = tau_block - D;

        % 每列对应一个输出点，与 v3 使用相同的 sinc 核及 L2 归一化。
        H = sinc(m - delta);
        H = H ./ vecnorm(H, 2, 1);

        % 始终从完整 x 中取样，允许读取相邻块的数据。
        Idx = n - D - m;
        valid_mask = (Idx >= 1) & (Idx <= N);
        X_mat = zeros(kernelLength, numel(n), 'like', x);
        X_mat(valid_mask) = x(Idx(valid_mask));

        y(n) = sum(X_mat .* H, 1).';

        % 释放当前块的大矩阵，避免下一块分配时仍保留上一块的数据。
        clear H Idx valid_mask X_mat
    end
end
