function  [RFT_paramenter_esti, RFT_no_paramenter_esti] = f_AFDM_RFT_v4(Rx_PS_remove_GPS, ...
    Tx_signal_PS_AFDM_eff, SNR_dB, true_delay, true_doppler, ...
    Min_delay_full, Max_delay_full, Min_doppler_full, ...
    Max_doppler_full, factor_delay_refine, factor_doppler_refine, cfg, para)
%% 版本信息
% 2026/09/16 V4.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：AFDM匹配滤波域做徙动补偿，然后实现相干累积。
%              1.按预设时延和多普勒网格做二维匹配滤波；
%              2.在二维匹配滤波域，用sinc函数插值实现分数时延、补偿基于多普勒计算出的分数距离徙动及慢时间相位；
%              3.做相干积累估计目标参数。
%
% Modification:
%   V3.0 - 实现完整 delay-Doppler 网格的时域相干 RFT；
%   V3.1 - 在函数内部根据 SNR、累计符号数和真实参数生成自适应搜索网格。
%   V4 - 修改RFT的实现逻辑，避免在每次迭代中都做一次二维匹配滤波，以降低RFT算法的计算复杂度。
%   V4.1 - 分离整数匹配滤波 delay 轴和分数 RFT delay 搜索轴；
%          在匹配滤波后的复数 delay 剖面上沿分数迁移轨迹插值。
%
% 输入参数：
%   Rx_echo_signal, Tx_signal_pulsedAFDM_eff - 接收回波和发射参考信号
%   para, cfg                               - 系统参数和配置参数
%   SNR_dB, true_delay, true_nor_doppler    - SNR 和真实目标参数
%   Min_delay_full, Max_delay_full           - 完整时延范围
%   Min_nor_Doppler_full, Max_nor_Doppler_full - 完整 Doppler 范围
%
% 输出参数：
%   RFT_map                   - delay-Doppler 相干积累功率面
%   est_delay                 - 估计归一化时延
%   est_nor_doppler           - 估计归一化 Doppler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 输入AFDM有效信号矩阵，已经由f_AFDM_Rx_preprocess_v1去掉前后保护区，每列只保留有效AFDM符号+tail，因此 RM 时间间隔也按该长度计算。
    samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;
    samples_per_symbol_up = samples_per_symbol * para.L;
    num_symbol = para.Nsym;
    % RFT 输出搜索网格：delay 可以为分数，Doppler 也可以为任意实数。SNR_dB
    [delay_search_grid, doppler_search_grid] = f_RFT_searchRange_v1(SNR_dB, true_delay, ...
        true_doppler, num_symbol, Min_delay_full, Max_delay_full, Min_doppler_full, ...
        Max_doppler_full, factor_delay_refine, factor_doppler_refine);
    % delay_search_grid = Min_delay_full:factor_delay_refine:Max_delay_full;
    % doppler_search_grid = Min_doppler_full:factor_doppler_refine:Max_doppler_full;

    % tmp_delay = 0.25;
    % tmp_doppler = 0.25;
    % true_delay = round(true_delay,3);
    % true_doppler = round(true_doppler,3);
    % delay_search_grid = (true_delay-tmp_delay):0.001:(true_delay+tmp_delay);
    % doppler_search_grid = (true_doppler-tmp_doppler):0.001:(true_doppler+tmp_doppler);

    len_delay_search_grid = length(delay_search_grid);
    len_doppler_search_grid = length(doppler_search_grid);


    % 二维匹配滤波只在唯一、连续的整数 delay 轴上计算。该坐标轴与上面的
    % 分数 RFT 输出轴相互独立，后续通过 sincInterpComplex 在任意分数
    % delay_search_grid + migration_delay 位置读取复数匹配结果。
    mf_delay_axis = ceil(Min_delay_full):floor(Max_delay_full);
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
        % Two_dim_MF_cube(:, :, n_sym) = f_DelayDopEst_ML_timeDomain_v10(factor_delay_refine, factor_doppler_refine, ...
        %     Tx_signal_PS_AFDM_eff(:, n_sym), Rx_PS_remove_GPS(:, n_sym), Min_delay_full, Max_delay_full, ...
        %     Min_nor_doppler_full, Max_nor_Doppler_full, para.L, para.N);
        Two_dim_MF_cube(:, :, n_sym) = f_DelayDopEst_ML_timeDomain_v12(Tx_signal_PS_AFDM_eff(:, n_sym), ...
            Rx_PS_remove_GPS(:, n_sym), mf_delay_axis, doppler_search_grid, para.L, para.N);
    end

    % AFDM符号起点
    symbol_sample_relative = (0:num_symbol - 1) * samples_per_symbol_up;
    sinc_half_len = 12;
    % RFT算法
    for doppler_idx = 1:len_doppler_search_grid
        nor_doppler = doppler_search_grid(doppler_idx);
        % nor_doppler = doppler_search_grid(18);
        % 计算当前采样点下的徙动，这里使用的是走-停模型
        migration_delay = nor_doppler * para.delta_f / para.fc .* symbol_sample_relative;

        % 匹配滤波函数内部只补偿了符号内从0开始的 doppler，补偿符号间慢时间相位。
        slow_time_phase = exp(-1i * 2 * pi * nor_doppler / (para.Nc * para.L) .* symbol_sample_relative);

        % 相干累积结果
        coherent_sum = zeros(len_delay_search_grid, 1, 'like', Two_dim_MF_cube);
        % 非相干累积结果
        no_coherent_sum = zeros(len_delay_search_grid, 1);
        for symbol_idx = 1:num_symbol
            % 当前 doppler 候选对应的复数 delay 剖面。
            delay_profile = Two_dim_MF_cube(:, doppler_idx, symbol_idx);

            % 正速度目标的峰值位于 initial_delay + migration_delay(symbol_idx)。
            % delay_search_grid 保持完整的分数值，小数部分由下面的 sinc插值直接处理，不需要预先拆成整数部分和小数部分。
            query_delay = delay_search_grid(:) + migration_delay(symbol_idx);

            % 在分数轨迹位置读取复数匹配滤波结果。
            trajectory_value = sincInterpComplex(mf_delay_axis, delay_profile, query_delay, sinc_half_len);

            coherent_sum = coherent_sum + trajectory_value * slow_time_phase(symbol_idx);
            no_coherent_sum = no_coherent_sum + abs(trajectory_value * slow_time_phase(symbol_idx)).^2;
        end
        % 相干累积
        RFT_complex(:, doppler_idx) = coherent_sum / num_symbol;
        % 非相干累积
        RFT_real(:, doppler_idx) = no_coherent_sum;
    end

    RFT_map = abs(RFT_complex).^2;
    RFT_map_no = RFT_real / num_symbol;

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
