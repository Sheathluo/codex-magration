function [set_nor_parameters_esti] = f_paraEstibyGradDesc_newtonOMP_single_trarget_v9(ML_value_refine_set_aver, ...
    AmbigFunc_pulsedAFDM_noDelayDop, Nor_delay_refineSearch_set, Nor_Doppler_refineSearch_set, N, L, ...
    L_delay, L_Doppler, max_norDelay, max_norDoppler, N_delay_peak, N_Doppler_peak, fs_up, delta_f, ...
    Num_iter_single, Num_iter_cyclic, Th_noTarget, Max_num_target, para, Migration_non, cfg)
%% 版本信息
% 2026/11/14 V9.0：仅支持单目标，先仿单目标性能，后续再改成多目标
% 2026/11/10 V9.0：配套加权模板和加权梯度v8。
% cfg.symbol_power_weight由主脚本按发射参考生成，原v5入口保持不变。
% 2026/07/17 V4.0
% 让这个函数能够匹配徙动目标。
% 2026/07/18 V5.0
% 在不改变 v4 计算流程的基础上，取消全部内部子函数，将上下文初始化、峰值检测、窗口生成、Newton refine
% 和目标响应重构全部展开到主函数中，便于逐行调试。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：基于 OMP/SIC 和 Newton 梯度 refine 的多目标归一化时延、
%              多普勒和功率估计。
% Modification：输入增加了不补偿的非相干累积结果，但是还是要使用补偿后的
%               累积结果作粗估计。修改时间：2026年07月17日；修改人：罗韬。
%
% 功能说明：
%   本函数输入粗距离徙动补偿后的二维匹配滤波功率面，迭代检测其中的目标峰值。
%   每检测到一个目标，先在峰值附近的小窗口内做 Newton refine，得到非整数
%   归一化时延、多普勒和目标功率估计；随后对所有已检测目标做循环 refinement，
%   并通过 SIC（Successive Interference Cancellation）从原始功率面中扣除已检测
%   目标的理论模糊函数响应，更新残差后继续检测下一个目标。
%
% 算法流程：
%   1. 根据输入门限和经验下限确定最终检测门限。
%   2. 在当前残差匹配滤波面中寻找最大峰值，若低于门限则停止。
%   3. 以峰值为中心裁剪局部时延/多普勒窗口，做单目标 Newton refine。
%   4. 对所有已检测目标做循环 refinement，降低目标间互扰。
%   5. 重构并扣除所有已检测目标的匹配滤波响应，更新残差并进入下一轮 OMP。
%
% 输入参数：
%   ML_value_refine_set_aver        - 接收数据得到的二维平均匹配滤波功率面
%   AmbigFunc_pulsedAFDM_noDelayDop - pulse-shaped AFDM 理论模糊函数
%   Nor_delay_refineSearch_set      - 归一化时延搜索集合
%   Nor_Doppler_refineSearch_set    - 归一化多普勒搜索集合
%   N, L                            - AFDM 子载波数和过采样倍数
%   L_delay, L_Doppler              - 接收端时延/多普勒搜索过采样倍数
%   max_norDelay, max_norDoppler    - 重构目标回波时使用的最大归一化时延/多普勒
%   N_delay_peak, N_Doppler_peak    - 峰值附近用于 Newton refine 的局部窗口半宽
%   fs_up, delta_f                  - 上采样后采样率和子载波间隔
%   Num_iter_single                 - 单目标 Newton refine 最大迭代控制参数
%   Num_iter_cyclic                 - 多目标循环 refinement 轮数
%   Th_noTarget                     - 无目标检测门限
%   Max_num_target                  - OMP 最多检测目标数
%   para                             - 系统参数结构体，使用其中的载频 para.fc
%   Migration_non                    - 未做距离徙动补偿的非相干累积结果
%
% 输出参数：
%   set_nor_parameters_esti - 3 x Num_detected_target 估计矩阵
%                             第 1 行：归一化时延
%                             第 2 行：归一化多普勒
%                             第 3 行：目标匹配滤波功率/增益估计
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% 1. 输入检查及反复使用参数的初始化
    % 默认返回空集；只有检测到目标后才裁剪为 3 x Num_detected_target。
    % set_nor_parameters_esti = [];
    [Len_x, Len_y] = size(ML_value_refine_set_aver);

    % 以下变量对应 v4 中 build_refine_context 的内容，现直接在主函数初始化。
    Ts = 1 / fs_up;

    % 离线理论模糊函数的采样倍率。
    L_delay_noDelayDop = 2;
    L_Doppler_noDelayDop = 2;
    L_delay_aver = 1 / (Nor_delay_refineSearch_set(2) - Nor_delay_refineSearch_set(1));
    L_Doppler_aver = 1 / (Nor_Doppler_refineSearch_set(2) - Nor_Doppler_refineSearch_set(1));
    delay_grid_step = Nor_delay_refineSearch_set(2) - Nor_delay_refineSearch_set(1);
    Doppler_grid_step = Nor_Doppler_refineSearch_set(2) - Nor_Doppler_refineSearch_set(1);
    Delta_tau = 1 / fs_up / L_delay_aver / (L_delay_noDelayDop / L_delay_aver);
    Delta_v = delta_f / L_Doppler_aver / (L_Doppler_noDelayDop / L_Doppler_aver);

    % 重构目标响应后，只截取接收搜索范围对应的二维功率面。
    echo_rows = 1:L_delay_noDelayDop:max_norDelay * L * L_delay_noDelayDop;
    echo_cols = N * L * L_Doppler_noDelayDop / 2 - max_norDoppler * L * L_Doppler_noDelayDop / 2 + 1 : ...
        N * L * L_Doppler_noDelayDop / 2 + max_norDoppler * L * L_Doppler_noDelayDop / 2;

    % 考虑强目标消除后的残留误差，沿用 v1 的经验门限下限以抑制假目标。
    max_ML_value = max(ML_value_refine_set_aver(:));
    Th_noTarget = max(Th_noTarget, 0.0044 * max_ML_value);

    % 预分配最多 Max_num_target 个目标，避免每检测一个目标都动态扩展矩阵。
    set_nor_parameters_esti = nan(3, Max_num_target);
    set_index_delay_set_fracEst = zeros(Max_num_target, 2 * N_delay_peak + 1);
    set_index_Doppler_set_fracEst = zeros(Max_num_target, 2 * N_Doppler_peak + 1);

    ML_value_residual = ML_value_refine_set_aver;
    Num_detected_target = 0;                        % 用于记录检测的是第几个目标

    % 计算徙动的系数，如果要修改徙动的计算规则，则需要修改这个参数
    samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;
    gamma = samples_per_symbol * L * delta_f / para.fc;

    % 正常 Newton 步不做限制；仅当原始步长超过一个粗网格时触发保护，
    % 并将该异常更新缩放到半个网格以内，避免跳出当前主峰吸引域。
    newton_step_guard_threshold = [delay_grid_step; Doppler_grid_step];
    max_guarded_newton_step = 0.5 * newton_step_guard_threshold;

    %% 2. OMP 检测、Newton 精估计和循环 SIC
    while Num_detected_target < Max_num_target
        %% 2.1 在当前残差面中检测最大峰值
        [value_max_refine, max_idx] = max(ML_value_residual(:));
        [I_row_max_refine, I_col_max_refine] = ind2sub(size(ML_value_residual), max_idx);

        % 为了低SNR仿真
        % if value_max_refine <= Th_noTarget
        %     break;
        % end

        %% 2.2 构造峰值附近的局部 refine 窗口
        % 当峰值靠近边界时，只平移局部窗口中心；粗估计仍使用原峰值索引，防止索引越界
        I_row_window_center = min(max(I_row_max_refine, N_delay_peak + 1), Len_x - N_delay_peak);
        I_col_window_center = min(max(I_col_max_refine, N_Doppler_peak + 1), Len_y - N_Doppler_peak);
        % 建立索引
        index_delay_set_fracEst = I_row_window_center - N_delay_peak:I_row_window_center + N_delay_peak;
        index_Doppler_set_fracEst = I_col_window_center - N_Doppler_peak:I_col_window_center + N_Doppler_peak;

        nor_delay_est = Nor_delay_refineSearch_set(I_row_max_refine);
        nor_Doppler_est = Nor_Doppler_refineSearch_set(I_col_max_refine);
        % 粗补偿徙动后的匹配滤波功率面
        D_func = ML_value_residual(index_delay_set_fracEst, index_Doppler_set_fracEst);
        % 未补偿徙动的功率面
        D_func_migration = Migration_non(index_delay_set_fracEst, index_Doppler_set_fracEst);

        % figure;
        % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
        % mesh(X, Y, 10 * log10(D_func / max(max(D_func))));
        % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
        % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
        % zlabel('Normalized average squared MF');
        % set(gca, 'FontSize', 10);
        % title('截取出的粗补偿徙动后实际匹配滤波');
        %
        % figure;
        % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
        % mesh(X, Y, 10 * log10(D_func_migration / max(max(D_func_migration))));
        % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
        % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
        % zlabel('Normalized average squared MF');
        % set(gca, 'FontSize', 10);
        % title('未补偿徙动的实际匹配滤波');

        %% 2.3 对新检测目标执行单目标 Newton refine
        nor_tau_list = Nor_delay_refineSearch_set(index_delay_set_fracEst) / ...
            L_delay_aver * L_delay_noDelayDop;
        nor_v_list = Nor_Doppler_refineSearch_set(index_Doppler_set_fracEst) / ...
            L_Doppler_aver * L_Doppler_noDelayDop;
        % 仍以补偿后功率面的离散峰值作为 Newton 初值。前级补偿量已经包含在
        % v7 理论面中，不再额外用平均残余徙动修正初值，避免采样点级徙动时
        % 因逐符号近似造成系统性时延偏差。
        tau_prev = nor_delay_est;
        doppler_prev = nor_Doppler_est;

        P_prev = value_max_refine;
        max_newton_iter = max(1, Num_iter_single - 1);

        % fprintf('\n粗估计结果:\n');
        % fprintf('粗估计时延= %.4f \n', tau_prev);
        % fprintf('粗估计多普勒= %.4f \n', doppler_prev);
        % fprintf('粗估计幅度= %.4f \n', P_prev);

        for iter = 1:max_newton_iter
            % 梯度函数返回局部模型相对 delay/Doppler 的一阶和二阶信息。
            [Ja_matrix, He_matrix] = f_grad2DAmbigFunc_dealyDoppler_noParfor_newton_v8( ...
                tau_prev, doppler_prev, AmbigFunc_pulsedAFDM_noDelayDop, D_func, nor_tau_list, ...
                nor_v_list, index_delay_set_fracEst, index_Doppler_set_fracEst, Ts, delta_f, ...
                Delta_tau, Delta_v, N, L, L_delay, L_Doppler, L_delay_noDelayDop, ...
                L_Doppler_noDelayDop, max_norDelay, max_norDoppler, P_prev, para, D_func, cfg);

            %%%%%%%%%%%%%%%%%%%% 检查海森矩阵是否正定 %%%%%%%%%%%%%%%%%%%%%%
            % hessian_type = f_isPositiveDefiniteHessian(-He_matrix)
            % eig_H = eig(-He_matrix)

            % 使用线性方程求解代替显式矩阵求逆，提高数值稳定性。
            %%%%%%%%%%%%%%%%%%%% 优化方法一：全步长牛顿法 %%%%%%%%%%%%%%%%%%%%%%
            del_delayDoppler = He_matrix \ Ja_matrix;

            %%%%%%%%%%%%%%%%%%%% 优化方法二：牛顿法-硬约束步长 %%%%%%%%%%%%%%%%%%%%%%
            % 使用一个函数来判断，如果补偿超过一个时延或者多普勒网格，直接缩放搜索步长，避免发散。
            del_delayDoppler = f_guardNewtonDelayDopplerStep(del_delayDoppler, ...
                newton_step_guard_threshold, max_guarded_newton_step);

            tau_new = tau_prev - del_delayDoppler(1);
            doppler_new = doppler_prev - del_delayDoppler(2);

            % 估计目标幅度
            % 根据估计出的时延和多普勒计算时延徙动
            if cfg.flag_delayMigration == 3
                gamma = 0;
            end

            tau_new_migration = tau_new + (0:para.Nsym-1).' .* gamma * doppler_new;

            nor_delay_esti_upSam = tau_new_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
            nor_Doppler_esti_upSam = doppler_new * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);

            % echo_unit = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor(AmbigFunc_pulsedAFDM_noDelayDop, ...
            %     para.Nc, para.L, L_delay_noDelayDop, L_delay_noDelayDop, ...
            %     1, nor_delay_esti_upSam, nor_Doppler_esti_upSam);

            echo_unit = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v5( ...
                AmbigFunc_pulsedAFDM_noDelayDop, 1, nor_delay_esti_upSam, ...
                nor_Doppler_esti_upSam, para, echo_rows, echo_cols, cfg);

            echo_local = echo_unit(index_delay_set_fracEst, index_Doppler_set_fracEst);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % figure
            % plot(Nor_delay_refineSearch_set, 10*log10(abs(ML_value_residual(:,I_col_max_refine))./max(max(abs(echo_unit)))), '-b*');
            % hold on;
            % plot(Nor_delay_refineSearch_set, 10*log10(abs(echo_unit(:,I_col_max_refine))./max(max(abs(echo_unit)))), '--ro');
            % xline(cfg.delay_vec, 'g--');
            % grid on; axis([-inf inf -30 2]);
            % set(gca,'FontSize',10);
            % xlabel('Normalized delay $\tau$', 'Interpreter', 'latex');
            % ylabel('Normalized average squared MF');
            % legend('Received','Reconstructed', 'True delay','Location','south');
            % annotation('textbox', [0.45 0.8 0.1 0.1], 'String', ['$\hat \tau = $', num2str(tau_new,'%10.5f')], ...
            %     'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 10,'Interpreter','latex');
            % title('时延切片对比');
            %
            % figure
            % plot(Nor_Doppler_refineSearch_set, 10*log10(abs(ML_value_residual(I_row_max_refine,:))./max(max(abs(echo_unit)))), '-b*');
            % hold on;
            % plot(Nor_Doppler_refineSearch_set, 10*log10(abs(echo_unit(I_row_max_refine,:))./max(max(abs(echo_unit)))), '--ro');
            % xline(cfg.doppler_vec, 'g--');
            % grid on; axis([-inf inf -30 2]);
            % set(gca,'FontSize',10);
            % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
            % ylabel('Normalized average squared MF');
            % legend('Received','Reconstructed', 'True doppler','Location','south');
            % annotation('textbox', [0.45 0.8 0.1 0.1], 'String', ['$\hat \nu = $', num2str(doppler_new,'%10.5f')], ...
            %     'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 10,'Interpreter','latex');
            % title('多普勒切片对比');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % P_fit = sum(D_func .* echo_local, "all") / sum(echo_local.^2, "all");
            %
            % R_unit  = D_func - echo_local;
            % R_shape = D_func - P_fit * echo_local;
            %
            % figure;
            % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
            % % mesh(X, Y, 10 * log10(R_unit / max(max(R_unit))));
            % mesh(X, Y, R_unit);
            % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
            % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
            % zlabel('Normalized average squared MF');
            % set(gca, 'FontSize', 10);
            % % axis([-inf inf -inf inf -40 2]);
            % % title('构造的理论模糊函数');
            %
            % figure;
            % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
            % mesh(X, Y, R_shape);
            % % mesh(X, Y, 10 * log10(R_shape / max(max(R_shape))));
            % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
            % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
            % zlabel('Normalized average squared MF');
            % set(gca, 'FontSize', 10);
            % % title('接收回波匹配滤波结果');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            % mesh(X, Y, 10 * log10(echo_local / max(max(echo_local))));
            % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
            % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
            % zlabel('Normalized average squared MF');
            % set(gca, 'FontSize', 10);
            % title('移位后理论模糊函数');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % save_file = 'ERR_result.mat';
            % var_name  = 'data_cube';
            %
            % new_matrix = err;   % 每次的二维矩阵
            %
            % if exist(save_file, 'file')
            %     S = load(save_file);
            %     data_cube = S.data_cube;
            %     data_cube(:, :, end+1) = new_matrix;
            % else
            %     data_cube = new_matrix;
            % end
            %
            % save(save_file, var_name);
            %
            % S = load('ERR_result.mat');
            % data_cube = S.data_cube;
            % junzhi = mean(data_cube,3);
            %
            %
            % figure;
            % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
            % mesh(X, Y, junzhi);
            % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
            % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
            % zlabel('Normalized average squared MF');
            % set(gca, 'FontSize', 10);
            % title('残差均值');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            P_new = sum(D_func .* echo_local, "all") ./ norm(echo_local, 'fro')^2;
            % P_new = 1;
            cost = norm(D_func - P_new * echo_local, 'fro')^2;

            del_iter = norm([tau_new, doppler_new, P_new] - [tau_prev, doppler_prev, P_prev]);
            tau_prev = tau_new;
            doppler_prev = doppler_new;
            P_prev = P_new;

            % if mod(iter, 2) == 1
            % if mod(iter, 1) == 0
            %     fprintf('\n第%d轮迭代:\n', iter);
            %     fprintf('真实归一化时延= %.5f \n', cfg.delay_vec(1));
            %     fprintf('估计归一化时延= %.5f \n', tau_prev);
            %     fprintf('归一化时延估计误差 = %.5f \n', cfg.delay_vec(1) - tau_prev);
            %
            %     fprintf('真实归一化多普勒= %.5f \n', cfg.doppler_vec(1));
            %     fprintf('估计归一化多普勒= %.5f \n', doppler_prev);
            %     fprintf('归一化多普勒估计误差 = %.5f \n', cfg.doppler_vec(1) - doppler_prev);
            %
            %     fprintf('估计幅度= %.5f \n', P_prev);
            %     fprintf('三个估计误差的范数 = %.13f \n', del_iter);
            %     fprintf('损失值 = %.6f \n', cost);
            % end

            if del_iter < 1e-13
                break;
            end
        end

        target_para_esti = [tau_prev; doppler_prev; P_prev];
        % 按照随机产生的内容保留规定位数的小数
        target_para_esti = round(target_para_esti, cfg.decimals);
        % 保存新目标参数及其窗口索引，供后续循环 精估计 重复使用。
        Num_detected_target = Num_detected_target + 1;
        set_nor_parameters_esti(:, Num_detected_target) = target_para_esti;
        set_index_delay_set_fracEst(Num_detected_target, :) = index_delay_set_fracEst;
        set_index_Doppler_set_fracEst(Num_detected_target, :) = index_Doppler_set_fracEst;
        %{
        %% 2.4 循环精细估计：扣除其他目标后重新估计当前目标
        for num_cyclic = 1:Num_iter_cyclic
            for n = Num_detected_target:-1:1
                % 扣除强目标后的目标回波匹配滤波功率面
                echo_otherTarget_esti = zeros(Len_x, Len_y);

                for k = 1:Num_detected_target
                    if k ~= n
                        % 展开 reconstruct_target_echo：重构第 k 个其他目标。
                        nor_delay_esti = set_nor_parameters_esti(1, k);
                        nor_Doppler_esti = set_nor_parameters_esti(2, k);
                        P_esti = set_nor_parameters_esti(3, k);

                        tau_new_migration = nor_delay_esti + (0:para.Nsym-1).' .* gamma * nor_Doppler_esti;
                        nor_delay_esti_upSam = tau_new_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
                        % nor_delay_esti_upSam = nor_delay_esti * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
                        nor_Doppler_esti_upSam = nor_Doppler_esti * L_Doppler_aver * ...
                            (L_Doppler_noDelayDop / L_Doppler_aver);

                        echo_target_esti = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v5( ...
                            AmbigFunc_pulsedAFDM_noDelayDop, P_esti, nor_delay_esti_upSam, nor_Doppler_esti_upSam, ...
                            para, echo_rows, echo_cols, cfg);

                        echo_otherTarget_esti = echo_otherTarget_esti + echo_target_esti;
                    end
                end

                % 当前目标观测面 = 原始功率面 - 其他已检测目标的理论响应。
                echo_Target_n_esti = ML_value_refine_set_aver - echo_otherTarget_esti;
                index_delay_n = set_index_delay_set_fracEst(n, :);
                index_Doppler_n = set_index_Doppler_set_fracEst(n, :);
                D_func = echo_Target_n_esti(index_delay_n, index_Doppler_n);

                %% 对第 n 个已检测目标展开单目标 Newton refine
                nor_tau_list = Nor_delay_refineSearch_set(index_delay_n) / ...
                    L_delay_aver * L_delay_noDelayDop;
                nor_v_list = Nor_Doppler_refineSearch_set(index_Doppler_n) / ...
                    L_Doppler_aver * L_Doppler_noDelayDop;

                tau_prev = set_nor_parameters_esti(1, n);
                doppler_prev = set_nor_parameters_esti(2, n);
                P_prev = set_nor_parameters_esti(3, n);
                max_newton_iter = max(1, Num_iter_single - 1);

                for iter = 1:max_newton_iter
                    % 梯度函数返回局部模型相对 delay/Doppler 的一阶和二阶信息。
                    [Ja_matrix, He_matrix] = f_grad2DAmbigFunc_dealyDoppler_noParfor_newton_v8( ...
                        tau_prev, doppler_prev, AmbigFunc_pulsedAFDM_noDelayDop, D_func, ...
                        nor_tau_list, nor_v_list, index_delay_set_fracEst, index_Doppler_set_fracEst, ...
                        Ts, delta_f, Delta_tau, Delta_v, N, L, L_delay, L_Doppler, ...
                        L_delay_noDelayDop, L_Doppler_noDelayDop, max_norDelay, ...
                        max_norDoppler, P_prev, para, D_func, cfg);

                    % 使用线性方程求解代替显式矩阵求逆，提高数值稳定性。
                    %%%%%%%%%%%%%%%%%%%% 优化方法一：正常牛顿法 %%%%%%%%%%%%%%%%%%%%%%
                    del_delayDoppler = He_matrix \ Ja_matrix;
                    % del_delayDoppler = f_guardNewtonDelayDopplerStep(del_delayDoppler, ...
                    %     newton_step_guard_threshold, max_guarded_newton_step);

                    %%%%%%%%%%%%%%%%%%%% 优化方法二：Levenberg-Marquardt (LM)​ 方法 %%%%%%%%%%%%%%%%%%%%%%
                    % lambda = 0.001; % 初始阻尼系数
                    % He_lm = He_matrix + lambda * diag(diag(He_matrix)); % 或者 eye(size(He_matrix))
                    % del_delayDoppler = He_lm \ Ja_matrix;

                    tau_new = tau_prev - del_delayDoppler(1);
                    doppler_new = doppler_prev - del_delayDoppler(2);

                    % 估计目标幅度
                    % 根据估计出的时延和多普勒计算时延徙动
                    tau_new_migration = tau_new + (0:para.Nsym-1).' .* gamma * doppler_new;

                    % 展开单位功率目标响应的重构与局部功率更新。
                    % nor_delay_esti_upSam = tau_new * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
                    nor_delay_esti_upSam = tau_new_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
                    nor_Doppler_esti_upSam = doppler_new * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);

                    echo_unit = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v5( ...
                        AmbigFunc_pulsedAFDM_noDelayDop, 1, nor_delay_esti_upSam, nor_Doppler_esti_upSam, ...
                        para, echo_rows, echo_cols, cfg);
                    echo_local = echo_unit(index_delay_set_fracEst, index_Doppler_set_fracEst);
                    P_new = sum(sum(D_func .* echo_local)) ./ norm(echo_local, 'fro')^2;

                    del_iter = norm([tau_new, doppler_new, P_new] - [tau_prev, doppler_prev, P_prev]);
                    tau_prev = tau_new;
                    doppler_prev = doppler_new;
                    P_prev = P_new;

                    if mod(iter, 1) == 0
                        fprintf('\n 循环估计过程 \n');
                        fprintf('\n第%d轮迭代:\n', iter);
                        fprintf('真实归一化时延= %.5f \n', cfg.delay_vec(1));
                        fprintf('估计归一化时延= %.5f \n', tau_prev);
                        fprintf('归一化时延估计误差 = %.5f \n', cfg.delay_vec(1) - tau_prev);

                        fprintf('真实归一化多普勒= %.5f \n', cfg.doppler_vec(1));
                        fprintf('估计归一化多普勒= %.5f \n', doppler_prev);
                        fprintf('归一化多普勒估计误差 = %.5f \n', cfg.doppler_vec(1) - doppler_prev);

                        fprintf('估计幅度= %.5f \n', P_prev);
                        fprintf('三个估计误差的范数 = %.13f \n', del_iter);
                        fprintf('损失值 = %.6f \n', cost);
                    end


                    if del_iter < 1e-13
                        break;
                    end
                end

                set_nor_parameters_esti(:, n) = [tau_prev; doppler_prev; P_prev];
            end
        end

        %% 2.5 SIC：重构所有已检测目标并更新残差面
        echo_detectTarget_esti = zeros(Len_x, Len_y);
        for k = 1:Num_detected_target
            % 展开 reconstruct_target_echo：重构第 k 个已检测目标。
            nor_delay_esti = set_nor_parameters_esti(1, k);
            nor_Doppler_esti = set_nor_parameters_esti(2, k);
            P_esti = set_nor_parameters_esti(3, k);

            tau_new_migration = nor_delay_esti + (0:para.Nsym-1).' .* gamma * nor_Doppler_esti;
            nor_delay_esti_upSam = tau_new_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
            % nor_delay_esti_upSam = nor_delay_esti * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
            nor_Doppler_esti_upSam = nor_Doppler_esti * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);
            echo_unit = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v5( ...
                AmbigFunc_pulsedAFDM_noDelayDop, P_esti, nor_delay_esti_upSam, nor_Doppler_esti_upSam, ...
                para, echo_rows, echo_cols, cfg);

            echo_target_esti = echo_unit;
            echo_detectTarget_esti = echo_detectTarget_esti + echo_target_esti;
        end
        ML_value_residual = ML_value_refine_set_aver - echo_detectTarget_esti;

        % 若残差中没有超过门限的峰值，则认为没有新的可检测目标。
        max_ML_value_residual = max(ML_value_residual(:));
        if max_ML_value_residual <= Th_noTarget
            break;
        end
    end
     %}
    %% 3. 整理输出
    % 裁剪掉预分配但未使用的列，恢复与 v4 一致的输出尺寸。
    if Num_detected_target == 0
        set_nor_parameters_esti = [];
    else
        set_nor_parameters_esti = set_nor_parameters_esti(:, 1:Num_detected_target);
    end

end
