function [RFT_paramenter_esti, RFT_no_paramenter_esti, RFT_map, RFT_map_no] = f_AFDM_RFT_v6(Rx_PS_remove_GPS, ...
    Tx_signal_PS_AFDM_eff, delay_search_grid, doppler_search_grid, cfg, para)
% V6：批量计算频域相关与分数徙动，复用所有时延候选的参考频谱。
% 按窗口内模板能量归一化，消除候选之间平均白噪声功率的差异。
% 前两个输出分别为相干、非相干估计 [tau*fs_up, fd/delta_f]。
% 第三、四输出为对应功率面，行对应时延，列对应多普勒。
% cfg 保留在输入接口中，与 f_AFDM_RFT_v5 的调用方式一致。

    num_symbol = para.Nsym;
    num_reference = size(Tx_signal_PS_AFDM_eff, 1);
    num_receive = size(Rx_PS_remove_GPS, 1);
    len_delay_search_grid = numel(delay_search_grid);
    len_doppler_search_grid = numel(doppler_search_grid);

    RFT_complex = zeros(len_delay_search_grid, len_doppler_search_grid);
    RFT_map_no = zeros(len_delay_search_grid, len_doppler_search_grid);

    % 符号起点间隔包含 CP 和前后保护区，单位为过采样后的采样点。
    symbol_sample_relative = (0:num_symbol - 1) * para.samples_per_symbol * para.L;
    migration_per_doppler = para.delta_f / para.fc * symbol_sample_relative;

    % 在原始参考后补零，为初始时延及正、负徙动留出空间。
    max_shift = ceil(max(abs(doppler_search_grid)) * max(abs(migration_per_doppler)));
    fft_length = 2^nextpow2(num_receive + num_reference + ceil(max(abs(delay_search_grid))) + max_shift);
    delay_frequency = ifftshift(-floor(fft_length/2):ceil(fft_length/2)-1).' / fft_length;
    reference_energy = sum(abs(Tx_signal_PS_AFDM_eff).^2, 1);
    fast_time = (0:num_reference - 1).' / (para.Nc * para.L);

    % 接收频谱和初始时延的逆变换基只计算一次，供全部候选复用。
    receive_spectrum = fft(Rx_PS_remove_GPS, fft_length, 1);
    delay_basis = exp(2j * pi * (delay_search_grid(:) * delay_frequency.')) / fft_length;

    % 窗口内模板能量也是时延的傅里叶级数。其系数由参考频谱的
    % 非循环自相关与接收矩形窗口决定，因而无需逐时延生成时域模板。
    energy_frequency = (0:fft_length - 1).' / fft_length;
    window_sum = zeros(fft_length, 1);
    window_sum(1) = num_receive;
    window_sum(2:end) = exp(1j * pi * energy_frequency(2:end) * (num_receive - 1)) ...
        .* sin(pi * energy_frequency(2:end) * num_receive) ./ sin(pi * energy_frequency(2:end));
    energy_basis = exp(-2j * pi * (delay_search_grid(:) * energy_frequency.'));

    % 每个候选多普勒只生成一次参考频谱，所有时延、符号一起计算。
    for doppler_idx = 1:len_doppler_search_grid
        nor_doppler = doppler_search_grid(doppler_idx);
        migration_delay = nor_doppler * migration_per_doppler;

        % 模板乘正号多普勒相位，相关时取共轭，与 V5 的符号内补偿约定一致。
        reference_doppler = Tx_signal_PS_AFDM_eff .* exp(2j * pi * nor_doppler * fast_time);
        reference_spectrum = fft(reference_doppler, fft_length, 1);

        % 先在频域相乘完成相关，再用相移读取 d0 + migration_delay。
        matched_spectrum = receive_spectrum .* conj(reference_spectrum);
        migration_phase = exp(2j * pi * (delay_frequency * migration_delay));
        correlation_value = delay_basis * (matched_spectrum .* migration_phase);

        % 频率从负到正排列后做非循环自相关，保留真实的频率差，
        % 避免把分数时延下相差一个采样频率的项错误合并。
        ordered_spectrum = fftshift(reference_spectrum, 1);
        spectrum_autocorrelation = ifft(abs(fft(ordered_spectrum, 2 * fft_length, 1)).^2, [], 1);
        energy_coefficient = spectrum_autocorrelation(1:fft_length, :) .* window_sum;
        % 正、负频率差互为共轭，取两倍实部；零频系数先减半以免重复。
        energy_coefficient(1, :) = energy_coefficient(1, :) / 2;
        energy_migration_phase = exp(-2j * pi * (energy_frequency * migration_delay));
        observed_energy = 2 * real(energy_basis * (energy_coefficient .* energy_migration_phase)) / fft_length^2;

        % 与逐候选模板实现相同的能量归一化，统一白噪声输出方差。
        % noise_gain_scale = sqrt(reference_energy ./ observed_energy);
        noise_gain_scale = 1;
        trajectory_value = correlation_value .* noise_gain_scale / (para.Nc * para.L);

        % 模板已补偿符号内多普勒，此处补偿符号间慢时间相位。
        slow_time_phase = exp(-2j * pi * nor_doppler / (para.Nc * para.L) .* symbol_sample_relative);
        RFT_complex(:, doppler_idx) = trajectory_value * slow_time_phase.' / num_symbol;
        RFT_map_no(:, doppler_idx) = mean(abs(trajectory_value).^2, 2);
    end

    RFT_map = abs(RFT_complex).^2;

    % 相干积累与非相干积累分别选峰，返回对应的归一化参数。
    [~, peak_idx] = max(RFT_map(:));
    [delay_idx, doppler_idx] = ind2sub(size(RFT_map), peak_idx);
    RFT_paramenter_esti = [delay_search_grid(delay_idx), doppler_search_grid(doppler_idx)];

    [~, peak_idx] = max(RFT_map_no(:));
    [delay_idx, doppler_idx] = ind2sub(size(RFT_map_no), peak_idx);
    RFT_no_paramenter_esti = [delay_search_grid(delay_idx), doppler_search_grid(doppler_idx)];
end
