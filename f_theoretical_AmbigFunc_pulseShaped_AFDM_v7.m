function E_mat = f_theoretical_AmbigFunc_pulseShaped_AFDM_v7(p, N, L, c1, mu4, tau_list, v_list, tau0, v0)
% v6功能：计算脉冲成形AFDM的理论模糊函数
% v7优化版：深度向量化处理，利用矩阵相乘替代了所有的 n, m, k 嵌套循环，极大提升运算速度。删除多余代码，便于阅读。
% 修改者：LuoTao 日期：2026年05月11日
% compute_EX2 计算公式 E{|X(tau,v)|^2}
% 输入:
%   p    - 长度为 NL 的实数脉冲成形向量
%   N    - 符号数
%   L    - 过采样倍数 L
%   c1   - 参数 c1
%   mu4  - 参数 mu4
%   tau  - 时延 tau
%   v    - 多普勒 v
% 输出:
%   E_X_tau_v_sq - 结果

    numTau = length(tau_list);
    numDoppler = length(v_list);

    NL = N * L;
    % ================= 1. 预计算固定索引与相位矩阵 =================
    n_vec_col = (0:N-1)';
    m_vec_row = 0:N-1;
    y_vec_row = 0:NL-1;

    % 预计算 Part 2 中独立于多普勒 v 的相位常数项
    Phase_c1 = exp(1j * 4 * pi * c1 * n_vec_col * m_vec_row);

    % 预计算用于生成多普勒频移矩阵 E_shift 的循环索引
    idx_E = mod(y_vec_row + n_vec_col * L, NL) + 1;

    % 预处理脉冲成形向量，用于基于频域相位的亚像素时延平移
    freq_vec = ((0:NL-1) - NL/2).';
    p_fft_centered = fftshift(fft(p));

    % 初始化输出矩阵
    E_mat = zeros(numTau, numDoppler);

    % ================= 2. 开始二维遍历计算 =================
    parfor ti = 1:numTau
        tau = tau_list(ti);

        % 利用频域相移实现高精度的时域循环移位
        % p_Tx：发射信号到接收端的回波信号参考
        p_Tx = real(ifft(fftshift(p_fft_centered .* exp(-1i*2*pi/NL * tau0 * freq_vec))));

        % p_Rx：接收端匹配滤波的参考信号
        p_Rx = real(ifft(fftshift(p_fft_centered .* exp(-1i*2*pi/NL * tau  * freq_vec))));

        % 预计算当前 tau 下的滑动窗乘积矩阵 W_mat (尺寸: NL x N)
        % 取代原先最内层 k 循环中的 p_Tx_tmp(k_mod) * conj(p_Rx_tmp(ktau_mod))
        W_mat = zeros(NL, N);
        for m = 0:N-1
            W_mat(:, m+1) = p_Tx .* circshift(conj(p_Rx), -m*L);
        end

        % 预分配并行循环的行结果，保证 parfor 内存安全
        E_row = zeros(1, numDoppler);

        % 在多普勒维度进行并行遍历
        for vi = 1:numDoppler
            if mod(vi, 500) == 0
                fprintf('###doppler row %d/%d done in delay row %d \n', vi, numDoppler, ti);
            end
            v = v_list(vi);
            v_diff = v - v0;

            % ---------------- 核心计算 ----------------
            % 生成基础多普勒相位向量
            E_v = exp(-1j * 2 * pi * v_diff * (0:NL-1)' / NL);

            % 根据预计算的循环索引，快速生成 E_shift 矩阵 (尺寸: N x NL)
            E_shift = E_v(idx_E);

            % 【提速核心】：用一次矩阵乘法完美替代原始的 k 循环累加 (尺寸: N x N)
            S_matrix = E_shift * W_mat;

            % 乘以对应的载波相位偏移 phase_n
            phase_n = exp(1j * 2 * pi * v_diff * n_vec_col / N);
            inner_sum_mat = S_matrix .* phase_n;

            % ---------------- Part 1 ----------------
            % 对应原逻辑中 m = 0 的一列向量
            phase_part1 = exp(-1j * 2 * pi * v_diff * n_vec_col / N);
            part1 = abs(phase_part1.' * inner_sum_mat(:, 1) )^2;

            % ---------------- Part 2 ----------------
            % 组合相位矩阵并跨维度求和
            Phase_mat2 = Phase_c1 .* phase_part1;
            sum_nm_all = sum(Phase_mat2 .* inner_sum_mat, 1);
            part2 = sum(abs(sum_nm_all).^2) / N * (mu4 - 2);

            % ---------------- Part 3 ----------------
            % 直接求所有元素的模的平方和 (F范数)
            part3 = norm(inner_sum_mat, 'fro')^2;

            % 汇总单个栅格点结果
            E_row(vi) = part1 + part2 + part3;
        end

        % 赋值给主输出矩阵
        E_mat(ti, :) = E_row;

        fprintf('Delay row %d/%d done\n', ti, numTau);
    end
end
