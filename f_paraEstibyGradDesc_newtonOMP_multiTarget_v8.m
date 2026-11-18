function [set_nor_parameters_esti] = f_paraEstibyGradDesc_newtonOMP_multiTarget_v8(ML_value_refine_set_aver, ...
    AmbigFunc_pulsedAFDM_noDelayDop, Nor_delay_refineSearch_set, Nor_Doppler_refineSearch_set, N, L, ...
    L_delay, L_Doppler, max_norDelay, max_norDoppler, N_delay_peak, N_Doppler_peak, fs_up, delta_f, ...
    Num_iter_single, Num_iter_cyclic, Th_noTarget, Max_num_target, para, Migration_non, cfg_trial)
%% 版本信息
% 2026/07/17 V4.0
% 让这个函数能够匹配徙动目标。
% 2026/07/18 V5.0
% 在不改变 v4 计算流程的基础上，取消全部内部子函数，将上下文初始化、峰值检测、窗口生成、Newton refine
% 和目标响应重构全部展开到主函数中，便于逐行调试。
% 2026/08/24 V7.0
% 根据归一化时延和多普勒网格限制牛顿迭代步长，添加Armijo线性搜索确保收敛到平稳点。把Armijo线性搜索单独写成一个函数
% 2026/08/29 V8.0
% 适配不徙动时候优化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function
% Description：基于NOMP精确搜索的单目标归一化时延、多普勒和功率估计。
% 修改时间：2026年08月29日；修改人：罗韬。
% Modification：1.适配不徙动时候的优化
%
% 功能说明：
%   本函数输入粗距离徙动补偿后的二维匹配滤波功率面，迭代检测其中的目标峰值。
%   每检测到一个目标，先在峰值附近的小窗口内做 Newton refine，得到非整数归一化时延、多普勒和目标功率估计；
%   随后对所有已检测目标做循环精确搜索，
%   并通过 SIC（Successive Interference Cancellation,连续干扰消除）从原始功率面中扣除已检测目标的理论模糊函数响应，
%   更新残差后继续检测下一个目标。
%
% 算法流程：
%   1. 根据输入门限和经验下限确定最终检测门限。
%   2. 在当前残差匹配滤波面中寻找最大峰值，若低于门限则停止。
%   3. 以峰值为中心裁剪局部时延/多普勒窗口，做单目标 Newton迭代精搜索。
%   4. 对所有已检测目标做循环精搜索，降低目标间互扰。
%   5. 重构并扣除所有已检测目标的匹配滤波响应，更新残差并进入下一轮NOMP。
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
%   Num_iter_single                 - 单目标 Newton 迭代精搜索的最大迭代控制参数
%   Num_iter_cyclic                 - 多目标循环迭代精搜索轮数
%   Th_noTarget                     - 无目标检测门限
%   Max_num_target                  - NOMP 最多检测目标数
%   para                            - 系统参数结构体，使用其中的载频 para.fc
%   Migration_non                   - 未做距离徙动补偿的非相干累积结果
%   cfg_trial                       - 参数配置结构体
%
% 输出参数：
%   set_nor_parameters_esti - 3 x Num_detected_target 估计矩阵
%                             第 1 行：归一化时延
%                             第 2 行：归一化多普勒
%                             第 3 行：目标匹配滤波功率/增益估计
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1. 输入检查及反复使用参数的初始化
    % 默认返回空集；只有检测到目标后才裁剪为 3 x Num_detected_target。
    set_nor_parameters_esti = [];
    if Max_num_target <= 0
        return;
    end

    [Len_x, Len_y] = size(ML_value_refine_set_aver);

    % 检查局部精搜索窗口是否小于匹配滤波功率面的尺寸，避免索引越界。
    if Len_x < 2 * N_delay_peak + 1 || Len_y < 2 * N_Doppler_peak + 1
        error('Peak refinement window exceeds ML_value_refine_set_aver size.');
    end

    Ts = para.Ts;       % 采样间隔

    % 离线理论模糊函数的过采样倍数
    L_delay_noDelayDop = para.L;
    L_Doppler_noDelayDop = para.L;
    %
    L_delay_aver = 1 / (Nor_delay_refineSearch_set(2) - Nor_delay_refineSearch_set(1));
    L_Doppler_aver = 1 / (Nor_Doppler_refineSearch_set(2) - Nor_Doppler_refineSearch_set(1));
    % 归一化时延和多普勒网格步长
    delay_grid_step = Nor_delay_refineSearch_set(2) - Nor_delay_refineSearch_set(1);
    Doppler_grid_step = Nor_Doppler_refineSearch_set(2) - Nor_Doppler_refineSearch_set(1);

    Delta_tau = 1 / fs_up / L_delay_aver / (L_delay_noDelayDop / L_delay_aver);
    Delta_v = delta_f / L_Doppler_aver / (L_Doppler_noDelayDop / L_Doppler_aver);

    % 重构目标响应后，只截取接收搜索范围对应的二维功率面。
    echo_rows = 1 : L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop;       % 重构目标相应的时延截取范围
    echo_cols = N*L*L_Doppler_noDelayDop/2 - max_norDoppler*L*L_Doppler_noDelayDop/2 + 1 : ...
        N*L*L_Doppler_noDelayDop/2 + max_norDoppler*L*L_Doppler_noDelayDop/2;   % 多普勒截取范围

    % 考虑强目标消除后的残留误差，用经验门限下限以抑制假目标。二维CFAR检测目标
    max_ML_value = max(ML_value_refine_set_aver(:));
    Th_noTarget = max(Th_noTarget, 0.0044 * max_ML_value);

    % 预分配最多 Max_num_target 个目标，避免每检测一个目标都动态扩展矩阵。
    set_nor_parameters_esti = nan(3, Max_num_target);
    set_index_delay_set_fracEst = zeros(Max_num_target, 2 * N_delay_peak + 1);
    set_index_Doppler_set_fracEst = zeros(Max_num_target, 2 * N_Doppler_peak + 1);

    ML_value_residual = ML_value_refine_set_aver;
    Num_detected_target = 0;

    % 计算徙动的系数，若修改徙动的计算规则，则需修改此参数
    samples_per_symbol = para.Nc + para.Ncp + 2 * para.N_preSuffix;

    if cfg_trial.flag_delayMigration == 3
        gamma = 0;
    else
        gamma = samples_per_symbol * L * delta_f / para.fc;
    end


    % 归一化时延和多普勒分别以各自粗搜索网格为异常方向判据；仅在方向超限时生效，将二维方向整体缩放到半个网格，
    % 随后仍由 Armijo 判断是否接受。
    newton_step_guard_threshold = [delay_grid_step; Doppler_grid_step];
    max_guarded_newton_step = 0.5 * newton_step_guard_threshold;

    % Armijo 回溯线搜索参数：alpha 从 1 开始，每次缩小为 beta 倍。
    armijo_c1 = 1e-4;
    armijo_beta = 0.5;
    armijo_max_backtrack = 12;
    armijo_min_step = 2^(-armijo_max_backtrack);

    % 停止条件不再混合不同量纲的功率、时延和多普勒。只有cost
    % 改善已经很小且归一化参数步长也很小时，才认为局部优化收敛。
    relative_loss_tolerance = 1e-6;
    normalized_step_tolerance = 1e-3;

     %% 2. OMP 检测、Newton 精估计和循环 SIC
    while Num_detected_target < Max_num_target
        %% 2.1 在当前残差面中检测最大峰值
        [value_max_refine, max_idx] = max(ML_value_residual(:));
        [I_row_max_refine, I_col_max_refine] = ind2sub(size(ML_value_residual), max_idx);

        % 这里取消对价目表的抑制，目的是为了仿真低SNR条件下的性能
        % if value_max_refine <= Th_noTarget
        %     break;
        % end

        %% 2.2 构造峰值附近的局部精索搜窗口
        % 当峰值靠近边界时，只平移局部窗口中心；粗估计仍使用原峰值索引
        I_row_window_center = min(max(I_row_max_refine, N_delay_peak + 1), Len_x - N_delay_peak);
        I_col_window_center = min(max(I_col_max_refine, N_Doppler_peak + 1), Len_y - N_Doppler_peak);

        index_delay_set_fracEst = I_row_window_center - N_delay_peak:I_row_window_center + N_delay_peak;
        index_Doppler_set_fracEst = I_col_window_center - N_Doppler_peak:I_col_window_center + N_Doppler_peak;

        nor_delay_est = Nor_delay_refineSearch_set(I_row_max_refine);
        nor_Doppler_est = Nor_Doppler_refineSearch_set(I_col_max_refine);
        % 截取出的接收回波匹配滤波结果
        D_func = ML_value_residual(index_delay_set_fracEst, index_Doppler_set_fracEst);
        % 截取出的无徙动的接收回波匹配滤波结果
        D_func_migration = Migration_non(index_delay_set_fracEst, index_Doppler_set_fracEst);

        figure;
        [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
        mesh(X, Y, 10 * log10(D_func / max(max(D_func))));
        ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
        xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
        zlabel('Normalized average squared MF');
        set(gca, 'FontSize', 10);
        title('补偿徙动后匹配滤波');

        figure;
        [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
        mesh(X, Y, 10 * log10(D_func_migration / max(max(D_func_migration))));
        ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
        xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
        zlabel('Normalized average squared MF');
        set(gca, 'FontSize', 10);
        title('存在徙动的匹配滤波');

        figure;
        [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
        mesh(X, Y, D_func);
        ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
        xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
        zlabel('Normalized average squared MF');
        set(gca, 'FontSize', 10);
        title('补偿徙动后匹配滤波');

        %% 2.3 对新检测目标执行单目标 Newton 精搜索
        nor_tau_list = Nor_delay_refineSearch_set(index_delay_set_fracEst) / L_delay_aver * L_delay_noDelayDop;
        nor_v_list = Nor_Doppler_refineSearch_set(index_Doppler_set_fracEst) / L_Doppler_aver * L_Doppler_noDelayDop;
        % 仍以补偿后功率面的离散峰值作为 Newton 迭代的初值。前级补偿量已经包含在理论面中，
        % 不再额外用平均残余徙动修正初值，避免采样点级徙动时因逐符号近似造成系统性时延偏差。
        tau_prev = nor_delay_est;
        doppler_prev = nor_Doppler_est;
        P_prev = value_max_refine;
        max_newton_iter = max(1, Num_iter_single - 1);

        fprintf('\n粗估计结果:\n');
        fprintf('粗估计时延= %.4f \n', tau_prev);
        fprintf('粗估计多普勒= %.4f \n', doppler_prev);
        fprintf('粗估计幅度= %.4f \n', P_prev);
        for iter = 1:max_newton_iter
            %% 在当前时延/多普勒处重构单位功率理论响应，并以最小二乘更新P。后根据当前时延/多普勒构造理论模糊函数与接收回波匹配滤波结果计算残差
            % 计算当前多普勒的徙动
            tau_last_migration = tau_prev + (0:para.Nsym-1).' .* gamma * doppler_prev;
            % 计算上采样后归一化时延/多普勒值
            nor_delay_esti_upSam = tau_last_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
            nor_Doppler_esti_upSam = doppler_prev * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);
            % 根据估计出的时延、多普勒及徙动构造理论模糊函数
            echo_unit_last = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v3( ...
                AmbigFunc_pulsedAFDM_noDelayDop, 1, nor_delay_esti_upSam, nor_Doppler_esti_upSam, ...
                para, echo_rows, echo_cols, cfg_trial);

            % 截取构造出的理论模糊函数的相应位置
            echo_local_last = echo_unit_last(index_delay_set_fracEst, index_Doppler_set_fracEst);
            % 计算接收回波与估计的时延/多普勒值构造出的理论模糊函数的残差
            cost_last = norm(D_func - P_prev * echo_local_last, 'fro')^2;

            % 梯度函数返回局部模型相对 delay/Doppler 的一阶和二阶信息。支持徙动与不徙动的优化
            if cfg_trial.flag_delayMigration == 3
                 [Ja_matrix, He_matrix] = f_grad2DAmbigFunc_dealyDoppler_noParfor_newton_v3( ...
                    tau_prev, doppler_prev, AmbigFunc_pulsedAFDM_noDelayDop, D_func, nor_tau_list, ...
                    nor_v_list, index_delay_set_fracEst, index_Doppler_set_fracEst, Ts, delta_f, ...
                    Delta_tau, Delta_v, N, L, L_delay, L_Doppler, L_delay_noDelayDop, L_Doppler_noDelayDop, ...
                    max_norDelay, max_norDoppler, P_prev);
            else
                [Ja_matrix, He_matrix] = f_grad2DAmbigFunc_dealyDoppler_noParfor_newton_v7( ...
                    tau_prev, doppler_prev, AmbigFunc_pulsedAFDM_noDelayDop, D_func, nor_tau_list, ...
                    nor_v_list, index_delay_set_fracEst, index_Doppler_set_fracEst, Ts, delta_f, ...
                    Delta_tau, Delta_v, N, L, L_delay, L_Doppler, L_delay_noDelayDop, ...
                    L_Doppler_noDelayDop, max_norDelay, max_norDoppler, P_prev, para, D_func, cfg_trial);
            end
            % 这里的参数向量为 x = [tau; Doppler]，Armijo 只搜索 x 的步长；
            % 功率 P 不沿 Newton 方向更新，而是在每个候选 x 处重新最小二乘估计。
            cost_gradient = -Ja_matrix;
            cost_hessian = -He_matrix;
            %%%%%%%%%%%%%%%%%%%% 优化方法一：全步长牛顿法 %%%%%%%%%%%%%%%%%%%%%%%
            % 判断海森矩阵是否正定
            cost_hessian_type = f_isPositiveDefiniteHessian(cost_hessian);
            if strcmp(cost_hessian_type, '正定')
                % 正定时用牛顿迭代方向 d = -Ja_matrix \ He_matrix;
                % 使用线性方程求解代替显式矩阵求逆，提高数值稳定性。
                parameter_direction = -(cost_hessian \ cost_gradient);
            elseif strcmp(cost_hessian_type, '半正定')
                % 正定时用牛顿迭代方向 d = -Ja_matrix \ He_matrix;
                parameter_direction = -(cost_hessian \ cost_gradient);
            else
                % 其他情况使用梯度方向
                parameter_direction = -cost_gradient;
            end

            % 方向导数 g_F^T*d 必须小于0，才说明沿 d 的小步更新会降低目标函数。
            % 数值误差或病态 Hessian 可能破坏这一条件，此时再次回退到负梯度。
            directional_derivative = real(cost_gradient.' * parameter_direction);
            % 判断优化方向和方向导数是不是病态，判断方向导数是否大于0。如果满足下列任意条件直接回退到负梯度方向更新
            if ~all(isfinite(parameter_direction)) || ~isfinite(directional_derivative) || directional_derivative >= 0
                parameter_direction = -cost_gradient;
                directional_derivative = real(cost_gradient.' * parameter_direction);
            end

            if ~all(isfinite(parameter_direction)) || ~isfinite(directional_derivative) || directional_derivative >= 0
                % 无可用下降方向时停止当前目标的 Newton refine，保持上一估计。
                break;
            end

            %%%%%%%%%%%%%%%%%%%% 优化方法二：牛顿法-硬约束步长 %%%%%%%%%%%%%%%%%%%%%%
            % 使用一个函数来判断，如果步长超过一个时延或者多普勒网格，直接缩放搜索步长，避免发散。
            parameter_direction = -(cost_hessian \ cost_gradient);
            parameter_direction = f_guardNewtonDelayDopplerStep(parameter_direction, ...
                newton_step_guard_threshold, max_guarded_newton_step);
            % 方向导数必须<0，才说明沿着该方向导数更新会降低cost
            directional_derivative = real(cost_gradient.' * parameter_direction);

            %%%%%%%%%%%%%%%%%%%% 优化方法三：Levenberg-Marquardt (LM)​ 方法 %%%%%%%%%%%%%%%%%%%%%%
            % lambda = 0.001; % 初始阻尼系数
            % He_lm = He_matrix + lambda * diag(diag(He_matrix)); % 或者 eye(size(He_matrix))
            % del_delayDoppler = He_lm \ Ja_matrix;

            %%%%%%%%%%%%%%%%%%%%%%%%%% 直接更新 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 直接更新 根据梯度计算归一化时延和多普勒 %%%% 下降方向应该是梯度的反方向
            % tau_new = tau_prev + parameter_direction(1);
            % doppler_new = doppler_prev + parameter_direction(2);
            % % 根据候选归一化时延和多普勒逐符号计算时延徙动，并在该候选参数处重构单位功率理论响应 A(tau_new, Doppler_new)。
            % tau_new_migration = tau_new + (0:para.Nsym-1).' .* gamma * doppler_new;
            % nor_delay_esti_upSam = tau_new_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
            % nor_Doppler_esti_upSam = doppler_new * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);
            %
            % echo_unit = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v3( ...
            %     AmbigFunc_pulsedAFDM_noDelayDop, 1, nor_delay_esti_upSam, ...
            %     nor_Doppler_esti_upSam, para, echo_rows, echo_cols, cfg_trial);
            % echo_local = echo_unit(index_delay_set_fracEst, index_Doppler_set_fracEst);
            %
            % % 固定候选时延和多普勒后，解析求解使残差最小的功率P_new，再计算剖面loss F_new = ||D_func-P_new*A_new||_F^2。
            % P_new = sum(D_func .* echo_local, 'all') / norm(echo_local, 'fro')^2;
            % cost_new = norm(D_func - P_new * echo_local, 'fro')^2;
            %%%%%%%%%%%%%%%%%%%%%%%%%% 直接更新 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%% Armijo线性搜索 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Armijo线性搜索。若候选点不满足充分下降条件，则按照alpha <- beta*alpha 回溯，直到接受或达到停止回溯条件。
            armijo_step = 1;
            armijo_accepted = false;
            for armijo_iter = 1:armijo_max_backtrack
                % 当前候选点 x_trial = x_current + alpha*d
                tau_new = tau_prev + armijo_step * parameter_direction(1);
                doppler_new = doppler_prev + armijo_step * parameter_direction(2);

                % 根据候选归一化时延和多普勒逐符号计算时延徙动，并在该候选参数处重构单位功率理论响应 A(tau_new, Doppler_new)。
                tau_new_migration = tau_new + (0:para.Nsym-1).' .* gamma * doppler_new;
                nor_delay_esti_upSam = tau_new_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
                nor_Doppler_esti_upSam = doppler_new * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);

                echo_unit = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v3( ...
                    AmbigFunc_pulsedAFDM_noDelayDop, 1, nor_delay_esti_upSam, ...
                    nor_Doppler_esti_upSam, para, echo_rows, echo_cols, cfg_trial);
                echo_local = echo_unit(index_delay_set_fracEst, index_Doppler_set_fracEst);

                figure;
                [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
                mesh(X, Y, echo_local);
                % mesh(X, Y, 10 * log10(echo_local / max(max(echo_local))));
                ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
                xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
                zlabel('Normalized average squared MF');
                set(gca, 'FontSize', 10);
                title('构造的理论模糊函数');

                % 固定候选时延和多普勒后，解析求解使残差最小的功率P_new，再计算剖面loss F_new = ||D_func-P_new*A_new||_F^2。
                P_new = sum(D_func .* echo_local, 'all') / norm(echo_local, 'fro')^2;
                cost_new = norm(D_func - P_new * echo_local, 'fro')^2;

                % Armijo 充分下降条件：F_new <= F_current + c1*alpha*(g_F^T*d)。由于 g_F^T*d<0，
                % 右侧小于当前 loss，因此被接受的候选点一定获得足够的下降。
                if isfinite(cost_new) && cost_new <= cost_last + armijo_c1 * armijo_step * directional_derivative
                    armijo_accepted = true;
                    break;
                end


                % 当前候选点下降不足，缩短步长后重新构造理论响应并计算 loss。
                armijo_step = armijo_beta * armijo_step;
                if armijo_step < armijo_min_step
                    fprintf('当前点下降不足，终止迭代 \n');
                    break;
                end
            end

            if ~armijo_accepted
                % 最大回溯次数内没有找到充分下降点：拒绝本轮参数更新，保留 tau_prev、doppler_prev 和 P_prev，避免继续迭代导致发散。
                fprintf('\n最大回溯次数内未找到充分下降点，拒绝第%d轮迭代的更新。\n', iter);
                fprintf('时延 = %.6f, 多普勒 = % .6f。\n', tau_new, doppler_new);
                break;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%% Armijo线性搜索 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % 计算此次迭代的改善
            relative_loss_improvement = (cost_last - cost_new) / max(cost_last, eps);
            normalized_parameter_step = norm([(tau_new - tau_prev) / delay_grid_step; ...
                (doppler_new - doppler_prev) / Doppler_grid_step]);
            del_iter = norm([tau_new, doppler_new, P_new] - [tau_prev, doppler_prev, P_prev]);

            % 更新参数
            tau_prev = tau_new;
            doppler_prev = doppler_new;
            P_prev = P_new;

            % if mod(iter, 2) == 1
            if mod(iter, 1) == 0
                fprintf('\n第%d轮迭代:\n', iter);

                fprintf('真实时延= %.5f \n', cfg_trial.delay_i_cal);
                fprintf('估计时延= %.5f \n', tau_prev);
                fprintf('时延估计误差 = %.5f \n', cfg_trial.delay_i_cal - tau_prev);

                fprintf('估计多普勒= %.5f \n', doppler_prev);
                fprintf('真实多普勒= %.5f \n', cfg_trial.doppler_i_cal);
                fprintf('多普勒估计误差 = %.5f \n', cfg_trial.doppler_i_cal - doppler_prev);

                fprintf('估计幅度= %.5f \n', P_prev);
                fprintf('三个估计误差的范数 = %.13f \n', del_iter);
                fprintf('损失值 = %.6f \n', cost_new);

            end

            if relative_loss_improvement < relative_loss_tolerance && ...
                    normalized_parameter_step < normalized_step_tolerance

                fprintf('相对loss改善及参数更新小于容限，终止迭代 \n');
                break;
            end

            % if del_iter < 1e-10
            %     break;
            % end

        end

        target_para_esti = [tau_prev; doppler_prev; P_prev];

        % 保存新目标参数及其窗口索引，供后续循环 精估计 重复使用。
        Num_detected_target = Num_detected_target + 1;
        set_nor_parameters_esti(:, Num_detected_target) = target_para_esti;
        set_index_delay_set_fracEst(Num_detected_target, :) = index_delay_set_fracEst;
        set_index_Doppler_set_fracEst(Num_detected_target, :) = index_Doppler_set_fracEst;
        %% 3. 整理输出
        % 裁剪掉预分配但未使用的列，恢复与 v4 一致的输出尺寸。
        if Num_detected_target == 0
            set_nor_parameters_esti = [];
        else
            set_nor_parameters_esti = set_nor_parameters_esti(:, 1:Num_detected_target);
        end
    end