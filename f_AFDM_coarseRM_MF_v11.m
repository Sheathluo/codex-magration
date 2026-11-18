function [ML_value_refine_set_aver, est_dop_RM_max, Migration_non, migration_compensate_bin] = f_AFDM_coarseRM_MF_v11( ...
    X_pulseShaped_remove_preSuffix, Tx_signal_pulsedAFDM_eff, para, Nor_delay_Search_set, Nor_doppler_Search_set, ...
    factor_delay_refine, factor_doppler_refine, Min_delay, Max_delay, Min_nor_doppler, ...
    Max_nor_doppler, cfg_trial)
%%  版本信息
%%2026/10/22 V11.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：先对所有符号进行二维时域匹配滤波，再根据粗多普勒假设进行距离徙动补偿（时延维循环移位），最后非相干积累或分段相干累积
%
% Modification:
%   f_AFDM_coarseRM_MF_v8:         - 增加返回变量，返回因补偿徙动的逐列移动量; 修改日期：2026年07月21日 修改人：LuoTao
%   f_AFDM_coarseRM_MF_v9:         - 支持不徙动时不需要搜索参数;   修改日期：2026年08月29日 修改人：LuoTao
%   f_AFDM_coarseRM_MF_v10:        - 实现分段相干累积;  修改日期：2026年09月11日 修改人：LuoTao
%   f_AFDM_coarseRM_MF_v11:        - 减少函数计算量：提前根据采样率计算徙动，若不同多普勒之间的徙动完全一致，则只需要搜索一个多普勒;
%                                    修改日期：2026年10月22日 修改人：LuoTao

% 输入参数：
%   X_pulseShaped_remove_preSuffix - 接收端有效 AFDM 符号矩阵，列为符号索引
%   Tx_signal_pulsedAFDM_eff       - 发射端有效 AFDM 符号的 pulse-shaped 匹配滤波参考
%   para                           - AFDM 与雷达参数结构体
%   Nor_Doppler_refineSearch_set   - 用于粗距离徙动补偿的归一化多普勒搜索集合
%   factor_delay_refine            - 归一化时延搜索步进
%   factor_Doppler_refine          - 归一化多普勒搜索步进
%   Min_delay, Max_delay           - 归一化时延搜索范围
%   Min_nor_Doppler, Max_nor_Doppler - 归一化多普勒搜索范围
%   cfg_trial                      - 参数配置结构体
%
% 输出参数：
%   ML_value_refine_set_aver        - 粗距离徙动补偿后的平均匹配滤波功率面行对应归一化时延搜索点，列对应归一化多普勒搜索点
%   est_dop_RM_max                  - 选出的粗 RM 补偿归一化多普勒值，不作为最终多普勒精估计
%   Migration_non                   - 不补偿徙动的功率面
%   migration_compensate_bin        - 徙动导致的位移量，第一行表示距离徙动的位移量，第二行表示速度徙动的位移量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 输入矩阵已经由 f_AFDM_Rx_preprocess_v1 去掉前后保护区，每列只保留 CP + 有效 AFDM 符号，因此 RM 时间间隔也按该长度计算。
    samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;

    % 预分配结果向量
    len_delay = length(Min_delay:factor_delay_refine:Max_delay);
    len_doppler = length(Min_nor_doppler:factor_doppler_refine:Max_nor_doppler);
    est_doppler_len = length(Nor_doppler_Search_set);
    est_dop_RM_max = 0;
    % 提前定义不同多普勒补偿后的功率面
    ML_value_refine_set_tensor = zeros(len_delay, len_doppler, est_doppler_len);
    ML_value_seg_coher_tensor = zeros(len_delay, len_doppler, est_doppler_len);
    % 因补偿徙动导致的后续的逐列位移量（包含距离徙动和速度徙动的位移量）
    range_migration_bin_tmp = zeros([est_doppler_len, para.Nsym]);
    migration_compensate_bin = zeros([2, para.Nsym]);

    % 保存每次补偿徙动后的结果
    migration_compensate_results = zeros([size(ML_value_refine_set_tensor),]);

    % ============ 第一步：对所有符号分别进行二维匹配滤波 ============
    ML_value_all = zeros(len_delay, len_doppler, para.Nsym);
    for n_sym = 1:para.Nsym
        ML_value_all(:, :, n_sym) = f_DelayDopEst_ML_timeDomain_v10(factor_delay_refine, factor_doppler_refine, ...
            Tx_signal_pulsedAFDM_eff(:, n_sym), X_pulseShaped_remove_preSuffix(:, n_sym), ...
            Min_delay, Max_delay, Min_nor_doppler, Max_nor_doppler, para.L, para.N);
    end
    % [~, max_idx_ML] = max(ML_value_all(:));
    % [~, ~, max_page_ML] = ind2sub(size(ML_value_all), max_idx_ML);
    Migration_non = sum(abs(ML_value_all).^2, 3) ./ para.Nsym;

    % figure;
    % [X, Y] = meshgrid(Nor_doppler_Search_set, Nor_delay_Search_set);
    % mesh(X, Y, 10 * log10(abs(ML_value_all(:,:,1)) / max(max(abs(ML_value_all(:,:,1))))));
    % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
    % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
    % zlabel('Normalized average squared MF');
    % set(gca, 'FontSize', 10);
    % % axis([-inf inf -inf inf -40 2]);
    % title('接收回波匹配滤波结果-粗估计函数');
    %
    % figure;
    % [X, Y] = meshgrid(Nor_doppler_Search_set, Nor_delay_Search_set);
    % mesh(X, Y, 10 * log10(Migration_non / max(max(Migration_non))));
    % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
    % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
    % zlabel('Normalized average squared MF');
    % set(gca, 'FontSize', 10);
    % % axis([-inf inf -inf inf -40 2]);
    % title('接收回波匹配滤波结果-粗估计函数');

    % 如果设置不徙动，直接返回匹配滤波结果
    if cfg_trial.flag_delayMigration == 3
        ML_value_refine_set_aver = Migration_non;
        return;
    end

    % 剔除了那些徙动量一样的多普勒，以减少计算量
    Nor_doppler_Search_real = [];
    last_migration_bin = 0;
    for dop_idx = 1:est_doppler_len
        % 根据当前粗多普勒假设计算每个符号的时延偏移量（整数采样点）
        velocity_search = Nor_doppler_Search_set(dop_idx) .* para.delta_f * para.c / (2*para.fc);
        %%%% 时间向量
        time_vec = (0:para.Nsym-1) * para.L * samples_per_symbol; %%%% !! 暂时使用符号的第一个点作为整个符号的时延
        %%%% 双程距离变化
        delta_range = velocity_search * 2 * time_vec;
        %%%% 预补偿徙动变化
        migration_delta_tau_norm = delta_range / para.c;
        % 徙动采样点转换为时延的真实徙动量
        migration_delta_tau_time = migration_delta_tau_norm / para.fs;
        % 时延的真实徙动量除以时间分辨单元，得到跨越距离门的数量
        migration_compensate = migration_delta_tau_time / para.delta_tau / factor_delay_refine;
        if Nor_doppler_Search_set(dop_idx) >= 0
            range_migration_bin = floor(migration_compensate);
        else
            range_migration_bin = ceil(migration_compensate);
        end
        % 计算徙动量与前一次徙动量相比有没有区别，如果徙动量完全一致，那么不需要搜索这一次
        save_flag = sum(last_migration_bin - range_migration_bin);
        if save_flag == 0
            last_migration_bin = range_migration_bin;
            continue
        else
            last_migration_bin = range_migration_bin;
            Nor_doppler_Search_real = [Nor_doppler_Search_real,Nor_doppler_Search_set(dop_idx)];
        end
    end

    % 为了补偿慢时间的多普勒
    filter_half_length = para.N_tapFiler * para.L / 2;
    data_start_offset = (para.Ncp + para.N_preSuffix) * para.L;
    % matlab从1开始的
    symbol_data_start = filter_half_length + data_start_offset + (0:para.Nsym - 1) * samples_per_symbol * para.L;
    symbol_data_start = symbol_data_start + 1;

    % ============ 第二步：遍历粗多普勒假设，进行距离徙动补偿和非相干积累 ============
    est_doppler_len_real = numel(Nor_doppler_Search_real);
    for dop_idx = 1:est_doppler_len_real
        % dop_idx = 31; % 补偿值为0
        % dop_idx = 40;
        % 根据当前粗多普勒假设计算每个符号的时延偏移量（整数采样点）
        velocity_search = Nor_doppler_Search_real(dop_idx) .* para.delta_f * para.c / (2*para.fc);
        %%%% 时间向量
        % time_vec = (0:para.Nsym-1) * para.L * samples_per_symbol + floor(para.L * samples_per_symbol / 2);
        time_vec = (0:para.Nsym-1) * para.L * samples_per_symbol; %%%% !! 暂时使用符号的第一个点作为整个符号的时延
        %%%% 双程距离变化
        delta_range = velocity_search * 2 * time_vec;
        %%%% 预补偿徙动变化
        migration_delta_tau_norm = delta_range / para.c;
        % 徙动采样点转换为时延的真实徙动量
        migration_delta_tau_time = migration_delta_tau_norm / para.fs;
        % 时延的真实徙动量除以时间分辨单元，得到跨越距离门的数量
        migration_compensate = migration_delta_tau_time / para.delta_tau / factor_delay_refine;

        if Nor_doppler_Search_real(dop_idx) >= 0
            range_migration_bin = floor(migration_compensate);
        else
            range_migration_bin = ceil(migration_compensate);
        end
        range_migration_bin_tmp(dop_idx,:) = range_migration_bin;
        % 对每个符号的二维匹配滤波结果沿时延维进行循环移位（对齐峰值）
        ML_delay_aligned = zeros(size(ML_value_all));
        % ML_slow_time_doppler_aligned = zeros(size(ML_value_all));

        for n_sym = 1:para.Nsym
            % 反向移动：目标时延应减去偏移量，故取负
            shift_amount = -range_migration_bin(n_sym);
            % 使用 circshift 沿第一维（时延）循环移位
            ML_delay_aligned(:, :, n_sym) = circshift(ML_value_all(:, :, n_sym), shift_amount, 1);
            % 若希望避免循环绕回，可用下面的补零平移
            % ML_aligned(:, :, n_sym) = zeroPadShift(ML_value_all(:, :, n_sym), shift_amount);
        end

        %%%%%%%%%%%%%%%%%%%%%%% 非相干累积 %%%%%%%%%%%%%%%%%%%%%%%%%%
        % 非相干积累：幅度平方后平均（！！！非相干累积就是会放大噪声，特别是在低信噪比条件下对目标幅度的估计性能会严重下降）
        ML_value_refine_set_aver_ = sum(abs(ML_delay_aligned).^2, 3) ./ para.Nsym;
        ML_value_refine_set_tensor(:, :, dop_idx) = ML_value_refine_set_aver_;

        %%%%%%%%%%%%%%%%%%%%%%% 相干累积 %%%%%%%%%%%%%%%%%%%%%%%%%%
        % 补偿慢时间的相位
        slow_time_phase = exp(-1i * 2 * pi * Nor_doppler_Search_real(dop_idx) / (para.Nc * para.L) .* symbol_data_start);
        ML_slow_time_doppler_aligned =  ML_delay_aligned .* reshape(slow_time_phase, 1, 1, para.Nsym);
         % 相干积累
        ML_value_coh = abs(sum(ML_slow_time_doppler_aligned,3)).^2./para.Nsym^2;
        % 相干累积符号数
        num_sym_coherent = 5;
        seg_count = ceil(para.Nsym/num_sym_coherent);
        % 提前生成分段累积后的存储位置
        [row, col] = size(ML_slow_time_doppler_aligned, 1, 2);
        seg_coherent_acc = zeros([row, col, seg_count]);
        for seg_idx = 1:seg_count
            if seg_idx == seg_count
                idx = num_sym_coherent*(seg_idx-1)+1:para.Nsym;
                seg_coherent_acc(:,:,seg_idx) = abs(sum(ML_slow_time_doppler_aligned(:,:,idx),3)).^2./(length(idx))^2;
            else
                idx = num_sym_coherent*(seg_idx-1)+1:num_sym_coherent*seg_idx;
                seg_coherent_acc(:,:,seg_idx) = abs(sum(ML_slow_time_doppler_aligned(:,:,idx),3)).^2./num_sym_coherent^2;
            end
        end
        ML_value_refine_seg = sum(abs(seg_coherent_acc).^2, 3) ./ seg_count;
        ML_value_seg_coher_tensor(:,:,dop_idx) = ML_value_refine_seg;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % figure;
        % [X, Y] = meshgrid(Nor_doppler_Search_set, Nor_delay_Search_set);
        % mesh(X, Y, 10 * log10(ML_value_refine_set_aver_ / max(max(ML_value_refine_set_aver_))));
        % % mesh(X, Y, ML_value_refine_set_aver);
        % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
        % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
        % zlabel('Normalized average squared MF');
        % set(gca, 'FontSize', 10);
        % % axis([-inf inf -inf inf -40 2]);
        % title('实际匹配滤波结果-非相干累积');
        %
        % figure;
        % mesh(X, Y, 10 * log10(ML_value_coh / max(max(ML_value_coh))));
        % % mesh(X, Y, ML_value_refine_set_aver_co);
        % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
        % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
        % zlabel('Normalized average squared MF');
        % set(gca, 'FontSize', 10);
        % axis([-inf inf -inf inf -60 2]);
        % title('实际匹配滤波结果-相干累积');
        %
        % figure;
        % mesh(X, Y, 10 * log10(ML_value_refine_seg / max(max(ML_value_refine_seg))));
        % % mesh(X, Y, ML_value_refine_set_aver_co);
        % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
        % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
        % zlabel('Normalized average squared MF');
        % set(gca, 'FontSize', 10);
        % % axis([-inf inf -inf inf -40 2]);
        % title('实际匹配滤波结果-分段相干累积');

    end

    %%%%%%%%%%%%%%%%%%%%%% 计算多普勒失配导致的功率损失 %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % M = para.Nsym;
    % P = samples_per_symbol * para.L;
    % D = para.Nc * para.L;
    %
    % delta_nu = 0.1;
    % residual_phase = exp(1i * 2*pi * delta_nu / D .* (0:M-1) * P);
    % coherent_power_gain = abs(sum(residual_phase) / M).^2;
    % loss_dB = 10*log10(coherent_power_gain);
    %
    % fprintf('相干功率损失 = %.2f dB\n', loss_dB);
    %%%%%%%%%%%%%%%%%%%%%% 计算多普勒失配导致的功率损失 %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ============ 第三步：按峰值锐度选择粗多普勒页面 ============
    % 选择全局最大值所在的粗补偿页，作为后续 NOMP/梯度精估计输入。
    [~, max_idx] = max(ML_value_refine_set_tensor(:));
    [~, ~, max_page] = ind2sub(size(ML_value_refine_set_tensor), max_idx);
    ML_value_refine_set_aver = ML_value_refine_set_tensor(:, :, max_page);
    est_dop_RM_max = Nor_doppler_Search_real(max_page);
    % 徙动导致的位移量，第一行表示距离徙动的位移量，第二行表示速度徙动的位移量
    migration_compensate_bin(1,:) = range_migration_bin_tmp(max_page,:);

    %{
    %%%%%%%%%%%%%%%%%%%%%%% 分段相干累积 %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 用非相干累积初步估计出时延和多普勒，然后再做分段相干累积输出
    [~, max_peak_idxs] = max(ML_value_refine_set_aver(:));
    [I_row_max_refine, I_col_max_refine] = ind2sub(size(ML_value_refine_set_aver), max_peak_idxs);
    nor_delay_est = Nor_delay_Search_set(I_row_max_refine);
    nor_doppler_est = Nor_doppler_Search_set(I_col_max_refine);

    %%%% 相干累积的基本原则：
    % 一、相干时间内目标移动距离不超过距离分辨率(速度)；
    % 二、相关时间内目标多普勒频移不超过多普勒分辨率(加速度)；
    % 三、相干时间目标的RCS起伏不超过阈值；
    % 根据距离分辨率计算相干累积的最长时间
    coherent_time_range = para.delta_range_rad/(abs(nor_doppler_est)*para.delta_f*para.lambda);
    % 根据多普勒分辨率计算相干累积的最长时间
    % coherent_time_doppler = para.delta_fd;
    time_AFDM_upsample = samples_per_symbol*para.L/para.fs_up;
    % num_sym_coherent = floor(coherent_time_range/time_AFDM_upsample);
    num_sym_coherent = 15;
    % 分段相干累积
    seg_count = ceil(para.Nsym/num_sym_coherent);
    % 提前生成分段累积后的存储位置
    [row, col] = size(ML_aligned, 1, 2);
    seg_coherent_acc = zeros([row, col, seg_count]);
    for seg_idx = 1:seg_count
        if seg_idx == seg_count
            idx = num_sym_coherent*(seg_idx-1)+1:para.Nsym;
            seg_coherent_acc(:,:,seg_idx) = abs(sum(ML_value_all(:,:,idx),3)).^2./(length(idx))^2;
        else
            idx = num_sym_coherent*(seg_idx-1)+1:num_sym_coherent*seg_idx;
            seg_coherent_acc(:,:,seg_idx) = abs(sum(ML_value_all(:,:,idx),3)).^2./num_sym_coherent^2;
        end
    end
    % 分段相干累积后再做相干累积
    ML_value_refine_set_aver_co = sum(abs(seg_coherent_acc).^2, 3) ./ seg_count;
    %%%%%%%%%%%%%%%%%%%%%%% 分段相干累积 %%%%%%%%%%%%%%%%%%%%%%%%%%
    %}

    % figure;
    % [X, Y] = meshgrid(Nor_doppler_Search_set, Nor_delay_Search_set);
    % mesh(X, Y, 10 * log10(ML_value_refine_set_aver / max(max(ML_value_refine_set_aver))));
    % % mesh(X, Y, ML_value_refine_set_aver);
    % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
    % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
    % zlabel('Normalized average squared MF');
    % set(gca, 'FontSize', 10);
    % title('最大值-补偿徙动后-非相干累积-粗估计函数');

    % figure;
    % [X, Y] = meshgrid(Nor_doppler_Search_set, Nor_delay_Search_set);
    % mesh(X, Y, 10 * log10(ML_value_refine_set_aver_co / max(max(ML_value_refine_set_aver_co))));
    % % mesh(X, Y, ML_value_refine_set_aver_co);
    % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
    % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
    % zlabel('Normalized average squared MF');
    % set(gca, 'FontSize', 10);
    % % axis([-inf inf -inf inf -80 2]);
    % title('分段相干累积-粗估计函数');


end

% ===== 辅助函数：补零平移（非循环） =====
function shifted = zeroPadShift(mat, shift)
    % mat: 二维矩阵（时延 × 多普勒）
    % shift: 整数，正数向下平移（时延增加），负数向上平移（时延减少）
    [N_delay, N_dopp] = size(mat);
    shifted = zeros(N_delay, N_dopp);
    if shift > 0
        % 向下平移：底部移出，顶部补零
        shifted(shift+1:end, :) = mat(1:end-shift, :);
    elseif shift < 0
        % 向上平移：顶部移出，底部补零
        shifted(1:end+shift, :) = mat(-shift+1:end, :);
    else
        shifted = mat;
    end
end
