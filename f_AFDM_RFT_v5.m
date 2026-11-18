function  [RFT_paramenter_esti, RFT_no_paramenter_esti] = f_AFDM_RFT_v5(Rx_PS_remove_GPS, ...
    Tx_signal_PS_AFDM_eff, delay_search_grid, doppler_search_grid, cfg, para)
% 逐符号做二维匹配滤波，再沿时延维用频域相移补偿距离徙动。
% 两个搜索网格由外部传入；分别输出相干、非相干估计 [tau*fs_up, fd/delta_f]。

    num_symbol = para.Nsym;

    len_delay_search_grid = length(delay_search_grid);
    len_doppler_search_grid = length(doppler_search_grid);

    % 匹配滤波用连续的整数时延轴，RFT 输出轴可以包含分数时延。
    mf_delay_axis = ceil(cfg.min_nor_delay):floor(cfg.max_nor_delay);
    len_mf_delay_axis = length(mf_delay_axis);

    % 预先分配向量
    % 二维匹配滤波结果
    Two_dim_MF_cube = zeros(len_mf_delay_axis, len_doppler_search_grid, num_symbol);
    % RFT校正距离徙动后结果
    % 相干累积结果
    RFT_complex = zeros(len_delay_search_grid, len_doppler_search_grid);
    % 非相干累积结果
    RFT_real = zeros(len_delay_search_grid, len_doppler_search_grid);

    % ============ 第一步：所有符号分别进行二维匹配滤波 ============
    for n_sym = 1:num_symbol
        Two_dim_MF_cube(:, :, n_sym) = f_DelayDopEst_ML_timeDomain_v12(Tx_signal_PS_AFDM_eff(:, n_sym), ...
            Rx_PS_remove_GPS(:, n_sym), mf_delay_axis, doppler_search_grid, para.L, para.N);
    end

    % 使用原始帧中相邻符号的起点间隔，包含 CP 和前后保护区。
    symbol_sample_relative = (0:num_symbol - 1) * para.samples_per_symbol * para.L;
    migration_per_doppler = para.delta_f / para.fc * symbol_sample_relative;

    % 跟 RFT_v2 一样，在时延轴末尾补零，为正、负徙动留出余量。
    max_shift = ceil(max(abs(doppler_search_grid)) * max(abs(migration_per_doppler)));
    fft_length = len_mf_delay_axis + 2 * max_shift;
    delay_frequency = ifftshift(-floor(fft_length/2):ceil(fft_length/2)-1).' / fft_length;

    % 在候选初始时延处直接计算逆 DFT，保留搜索网格的小数部分。
    delay_offset = delay_search_grid(:) - mf_delay_axis(1);
    idft_basis = exp(2j * pi * (delay_offset * delay_frequency.')) / fft_length;

    % ============ 第二步：按候选多普勒补偿徙动并积累 ============
    for doppler_idx = 1:len_doppler_search_grid
        nor_doppler = doppler_search_grid(doppler_idx);
        migration_delay = nor_doppler * migration_per_doppler;

        % 每列是一个符号的复数时延剖面；正号相移读取 d0 + migration_delay。
        delay_profiles = reshape(Two_dim_MF_cube(:, doppler_idx, :), len_mf_delay_axis, num_symbol);
        delay_spectrum = fft(delay_profiles, fft_length, 1);
        migration_phase = exp(2j * pi * (delay_frequency * migration_delay));
        trajectory_value = idft_basis * (delay_spectrum .* migration_phase);

        % 超出已计算的匹配滤波时延轴时，仍按原来的规则置零。
        query_delay = delay_search_grid(:) + migration_delay;
        trajectory_value(query_delay < mf_delay_axis(1) | query_delay > mf_delay_axis(end)) = 0;

        % 匹配滤波函数内部只补偿了符号内从0开始的 doppler，补偿符号间慢时间相位。
        slow_time_phase = exp(-1i * 2 * pi * nor_doppler / (para.Nc * para.L) .* symbol_sample_relative);

        RFT_complex(:, doppler_idx) = trajectory_value * slow_time_phase.' / num_symbol;
        RFT_real(:, doppler_idx) = mean(abs(trajectory_value).^2, 2);
    end

    RFT_map = abs(RFT_complex).^2;
    RFT_map_no = RFT_real;

    % figure;
    % [X, Y] = meshgrid(doppler_search_grid, delay_search_grid);
    % mesh(X, Y, 10 * log10(RFT_map_no / max(max(RFT_map_no))));
    % % mesh(X, Y, ML_value_refine_set_aver);
    % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
    % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
    % zlabel('Normalized average squared MF');
    % set(gca, 'FontSize', 10);
    % axis([-inf inf -inf inf -40 2]);
    % title('实际匹配滤波结果-相干累积');

    % 相干累积的参数估计
    [~, peak_idx] = max(RFT_map(:));
    [delay_idx, doppler_idx] = ind2sub(size(RFT_map), peak_idx);
    est_delay = delay_search_grid(delay_idx);
    est_doppler = doppler_search_grid(doppler_idx);
    RFT_paramenter_esti = [est_delay, est_doppler];

    % 非相干累积的参数估计
    [~, peak_idx] = max(RFT_map_no(:));
    [delay_idx, doppler_idx] = ind2sub(size(RFT_map_no), peak_idx);
    est_delay = delay_search_grid(delay_idx);
    est_doppler = doppler_search_grid(doppler_idx);
    RFT_no_paramenter_esti = [est_delay, est_doppler];

end
