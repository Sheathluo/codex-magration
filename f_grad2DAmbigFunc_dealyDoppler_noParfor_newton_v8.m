function [J, H, response_migration] = f_grad2DAmbigFunc_dealyDoppler_noParfor_newton_v8( ...
    tau0, v0, AFDM_AF_time, Rx_echo_MF, nor_tau_list, nor_v_list, index_delay_set_fracEst, ...
    index_Doppler_set_fracEst, Ts, delta_f, Delta_tau, Delta_v, N, L, ...
    L_delay, L_Doppler, L_delay_noDelayDop, L_Doppler_noDelayDop, ...
    max_norDelay, max_norDoppler, P_sig_strong_nor, para, D_func, cfg)
%% 版本信息
% 2026/07/22 V7.0
% 2026/11/10 V8.0：基于v7，对响应及全部一、二阶导数同步加入符号功率权重。
% 配套响应函数：f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v5。
% cfg.symbol_power_weight为发射参考能量平方权重，平均值1，迭代时固定。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function：计算经过逐符号二维徙动补偿后的理论模糊函数及其解析导数，
%           并构造归一化时延、归一化 Doppler 的 Newton 梯度和 Hessian。
%
% 本版本在 v6 的距离徙动理论模型上增加前级补偿算子：
%   1. 每个符号先生成 tau_n = tau0 + n*gamma*v0 对应的理论响应；
%   2. 将理论响应映射到与实际匹配滤波面相同的接收搜索网格；
%   3. 按 cfg_trial.migration_compensate_bin 对时延维和 Doppler 维移位；
%   4. 对移位后的理论面及其全部一、二阶导数做跨符号平均。
%
% cfg_trial.migration_compensate_bin 的定义：
%   第 1 行 - 前级记录的距离徙动 bin，v7 取其负值作为时延维实际移位量；
%   第 2 行 - 前级记录的 Doppler 徙动 bin，v7 取其负值作为 Doppler 维
%            实际移位量；当前尚未补偿 Doppler 徙动时该行应为全零。
%   每一列对应一个 AFDM 符号，单位分别是接收 ML 面的行、列搜索栅格。
%
% 注意：补偿位移由前级粗估计固定，在一次 Newton refinement 中不随
% tau0/v0 连续变化。因此模型和所有导数面执行相同的 circshift，不对
% floor/ceil 得到的离散位移本身求导。
%
% 输入参数：
%   tau0, v0                         - 当前归一化时延和归一化 Doppler
%   AFDM_AF_time                     - 完整的单位功率二维理论模糊函数
%   Rx_echo_MF                       - 兼容旧接口，本版本不直接使用
%   nor_tau_list, nor_v_list         - 兼容旧接口，本版本不直接使用
%   index_delay_set_fracEst          - Newton 局部窗口在接收时延轴上的索引
%   index_Doppler_set_fracEst        - Newton 局部窗口在接收 Doppler 轴上的索引
%   Ts                               - 上采样采样周期，单位 s
%   delta_f                          - 子载波间隔，单位 Hz
%   Delta_tau, Delta_v               - 理论模糊函数时延/Doppler 网格间隔
%   N, L                             - AFDM 子载波数和时域过采样倍数
%   L_delay, L_Doppler               - 接收搜索面的时延/Doppler 采样倍率
%   L_delay_noDelayDop               - 理论模糊函数时延维采样倍率
%   L_Doppler_noDelayDop             - 理论模糊函数 Doppler 维采样倍率
%   max_norDelay, max_norDoppler     - 接收搜索区域边界
%   P_sig_strong_nor                 - 当前目标匹配滤波功率/增益估计
%   para                             - 波形参数，符号数使用 para.Nsym
%   D_func                           - 当前 Newton 局部窗口的实际补偿后 ML 面
%   cfg                              - 当前配置，包含补偿位移和symbol_power_weight
%
% 输出参数：
%   J                                - 2 x 1 Newton 梯度向量
%   H                                - 2 x 2 Newton Hessian 矩阵
%   response_migration               - 可选输出；当前局部窗口内、经徙动补偿后的单位功率理论响应面 A。
%                                      用于数值微分诊断，不影响原有双输出调用。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% 1. 输入检查及补偿位移读取
    [num_delay, num_doppler] = size(AFDM_AF_time);
    symbol_weights = cfg.symbol_power_weight(1:para.Nsym);
    symbol_weights = symbol_weights(:);
    expected_delay = N * L * L_delay_noDelayDop;
    expected_doppler = N * L * L_Doppler_noDelayDop;

    % v8 粗补偿函数保存的是计算出的正向徙动 bin；实际对 ML 面采用反向移位。
    delay_shift_amount = -cfg.migration_compensate_bin(1, :).';
    doppler_shift_amount = -cfg.migration_compensate_bin(2, :).';
    compensation_shift = [delay_shift_amount, doppler_shift_amount];

    compensation_shift = round(compensation_shift);

    delay_bins = [0:(num_delay / 2 - 1), -num_delay / 2:-1];
    doppler_bins = [0:(num_doppler / 2 - 1), -num_doppler / 2:-1];

    %% 2. 构造 FFT 相移、导数系数和接收搜索区域
    doppler_frequency_spacing = Delta_v * L_Doppler_noDelayDop;
    delay_scale = -1i * 2 * pi / num_delay * Ts / Delta_tau;
    doppler_scale = -1i * 2 * pi / num_doppler * doppler_frequency_spacing / Delta_v;

    samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;
    if cfg.flag_delayMigration == 3
        gamma = 0;
    else
        gamma = samples_per_symbol * L * delta_f / para.fc;
    end

    symbol_slope = (0:para.Nsym-1).' * gamma;

    delay_derivative = (delay_scale .* delay_bins).';
    doppler_derivative = doppler_scale .* doppler_bins;
    doppler_phase = exp(doppler_scale * v0 .* doppler_bins);
    doppler_phase_d1 = doppler_derivative .* doppler_phase;
    doppler_phase_d2 = doppler_derivative.^2 .* doppler_phase;
    initial_delay_phase = exp(delay_scale * tau0 .* delay_bins).';
    AFDM_AF_freq = fft(AFDM_AF_time, [], 1);

    delay_stride = L_delay_noDelayDop / L_delay;

    delay_range = 1:delay_stride:(max_norDelay * L * L_delay_noDelayDop);
    doppler_center = num_doppler / 2;
    doppler_half_width = max_norDoppler * L * L_Doppler_noDelayDop / 2;
    doppler_range = (doppler_center - doppler_half_width + 1): ...
        (doppler_center + doppler_half_width);
    search_size = [numel(delay_range), numel(doppler_range)];

    %% 3. 按相同二维补偿位移分组，生成并移动理论面及全部导数面
    % 同一组内的符号具有相同 circshift。利用 FFT 和移位的线性性质，
    % 先累积组内相位再统一移位，与逐个符号生成 response_n 后移位完全等价。
    response_acc = zeros([search_size, 6]);
    unique_shift = unique(compensation_shift, 'rows', 'stable');

    for idx_shift = 1:size(unique_shift, 1)
        shift_curr = unique_shift(idx_shift, :);
        group_mask = compensation_shift(:, 1) == shift_curr(1) & ...
            compensation_shift(:, 2) == shift_curr(2);
        slope_group = symbol_slope(group_mask);
        weight_group = symbol_weights(group_mask);

        migration_phase_group = exp(slope_group * (delay_scale * v0 .* delay_bins));
        % w_k与候选时延/多普勒无关，三个相位和使用完全相同的权重。
        weighted_migration_phase_group = weight_group .* migration_phase_group;
        phase_sum = initial_delay_phase .* sum(weighted_migration_phase_group, 1).';
        slope_phase_sum = initial_delay_phase .* sum(slope_group .* weighted_migration_phase_group, 1).';
        slope2_phase_sum = initial_delay_phase .* sum(slope_group.^2 .* weighted_migration_phase_group, 1).';

        % 六页频域输入用于构造 A、A_tau、A_v、A_tau,tau、A_tau,v、A_v,v。
        delay_frequency_stack = cat(3, ...
            phase_sum .* AFDM_AF_freq, ...
            delay_derivative .* phase_sum .* AFDM_AF_freq, ...
            delay_derivative.^2 .* phase_sum .* AFDM_AF_freq, ...
            delay_derivative .* slope_phase_sum .* AFDM_AF_freq, ...
            delay_derivative.^2 .* slope_phase_sum .* AFDM_AF_freq, ...
            delay_derivative.^2 .* slope2_phase_sum .* AFDM_AF_freq);
        doppler_frequency_stack = fft(ifft(delay_frequency_stack, [], 1), [], 2);

        response_frequency_stack = cat(3, ...
            doppler_frequency_stack(:, :, 1) .* doppler_phase, ...
            doppler_frequency_stack(:, :, 2) .* doppler_phase, ...
            doppler_frequency_stack(:, :, 4) .* doppler_phase + ...
                doppler_frequency_stack(:, :, 1) .* doppler_phase_d1, ...
            doppler_frequency_stack(:, :, 3) .* doppler_phase, ...
            doppler_frequency_stack(:, :, 5) .* doppler_phase + ...
                doppler_frequency_stack(:, :, 2) .* doppler_phase_d1, ...
            doppler_frequency_stack(:, :, 6) .* doppler_phase + ...
                2 * doppler_frequency_stack(:, :, 4) .* doppler_phase_d1 + ...
                doppler_frequency_stack(:, :, 1) .* doppler_phase_d2);
        response_stack = real(ifft(response_frequency_stack, [], 2));

        % 位移量单位属于接收 ML 搜索网格，因此必须先抽取搜索区域再移位。
        response_search_stack = response_stack(delay_range, doppler_range, :);
        response_search_stack = circshift(response_search_stack, [shift_curr(1), shift_curr(2), 0]);
        response_acc = response_acc + response_search_stack;
    end

    response_acc = response_acc / para.Nsym;
    response_local = response_acc(index_delay_set_fracEst, index_Doppler_set_fracEst, :);

    response_migration = response_local(:, :, 1);
    response_first_delay = response_local(:, :, 2);
    response_first_doppler = response_local(:, :, 3);
    response_second_delay = response_local(:, :, 4);
    response_second_delay_doppler = response_local(:, :, 5);
    response_second_doppler = response_local(:, :, 6);


    % figure;
    % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
    % mesh(X, Y, 10 * log10(D_func / max(max(D_func))));
    % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
    % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
    % zlabel('Normalized average squared MF');
    % set(gca, 'FontSize', 10);
    % title('移位后实际匹配滤波');
    %
    % figure;
    % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
    % mesh(X, Y, 10 * log10(response_migration / max(max(response_migration))));
    % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
    % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
    % zlabel('Normalized average squared MF');
    % set(gca, 'FontSize', 10);
    % title('移位后理论模糊函数');
    %
    % err = 10 * log10(D_func / max(max(D_func))) - 10 * log10(response_migration / max(max(response_migration)));
    % max(max(abs(err)))
    % figure;
    % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
    % mesh(X, Y, err);
    % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
    % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
    % zlabel('Normalized average squared MF');
    % set(gca, 'FontSize', 10);
    % title('残差');

    %% 4. 使用补偿后的同一理论面及其导数构造 Newton 梯度和 Hessian
    % 最大化目标：Q = 2*P*<Y,A> - P^2*<A,A>。
    residual = D_func - P_sig_strong_nor * response_migration;

    gradient_delay = 2 * sum(residual .* response_first_delay * P_sig_strong_nor, 'all');
    gradient_doppler = 2 * sum(residual .* response_first_doppler * P_sig_strong_nor, 'all');

    hessian_delay_delay = -2 * ( ...
        sum(P_sig_strong_nor^2 * response_first_delay .* ...
        response_first_delay, 'all') - ...
        sum(residual .* response_second_delay * P_sig_strong_nor, 'all'));
    hessian_doppler_doppler = -2 * ( ...
        sum(P_sig_strong_nor^2 * response_first_doppler .* ...
        response_first_doppler, 'all') - ...
        sum(residual .* response_second_doppler * P_sig_strong_nor, 'all'));
    hessian_delay_doppler = -2 * ( ...
        sum(P_sig_strong_nor^2 * response_first_delay .* response_first_doppler, 'all') - ...
        sum(residual .* response_second_delay_doppler * P_sig_strong_nor, 'all'));

    J = [gradient_delay; gradient_doppler];
    H = [hessian_delay_delay, hessian_delay_doppler; ...
         hessian_delay_doppler, hessian_doppler_doppler];
end
