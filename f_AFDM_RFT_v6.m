function [RFT_paramenter_esti, RFT_no_paramenter_esti, RFT_map, RFT_map_no] = f_AFDM_RFT_v6(Rx_PS_remove_GPS, ...
    Tx_signal_PS_AFDM_eff, delay_search_grid, doppler_search_grid, cfg, para)
% V6：对发射参考做频域分数时延移位，在完整接收窗口内匹配并积累。
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

    % 逐多普勒候选构造模板，再沿候选运动轨迹匹配。
    for doppler_idx = 1:len_doppler_search_grid
        nor_doppler = doppler_search_grid(doppler_idx);
        migration_delay = nor_doppler * migration_per_doppler;

        % 模板乘正号多普勒相位，相关时取共轭，与 V5 的符号内补偿约定一致。
        reference_doppler = Tx_signal_PS_AFDM_eff .* exp(2j * pi * nor_doppler * fast_time);
        reference_spectrum = fft(reference_doppler, fft_length, 1);
        trajectory_value = zeros(len_delay_search_grid, num_symbol);

        for symbol_idx = 1:num_symbol
            query_delay = delay_search_grid(:).' + migration_delay(symbol_idx);

            % 负号频域相位将参考移至 d0 + migration_delay，保留分数采样点。
            delay_phase = exp(-2j * pi * (delay_frequency * query_delay));
            delayed_reference = ifft(reference_spectrum(:, symbol_idx) .* delay_phase, [], 1);
            observed_reference = delayed_reference(1:num_receive, :);

            % 每个候选使用接收窗口内的实际模板能量，统一白噪声输出方差。
            % 整数时延且模板完整落在窗口内时，等于 V5 的 1/(Nc*L) 相关值。
            observed_energy = sum(abs(observed_reference).^2, 1);
            noise_gain_scale = sqrt(reference_energy(symbol_idx) ./ observed_energy);
            correlation_value = Rx_PS_remove_GPS(:, symbol_idx).' * conj(observed_reference);
            trajectory_value(:, symbol_idx) = (correlation_value .* noise_gain_scale).' / (para.Nc * para.L);
        end

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
