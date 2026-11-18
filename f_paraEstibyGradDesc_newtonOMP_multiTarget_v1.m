function [set_nor_parameters_esti] = f_paraEstibyGradDesc_newtonOMP_multiTarget_v1(ML_value_refine_set_aver, ...
    AmbigFunc_pulsedAFDM_noDelayDop,Nor_delay_refineSearch_set,Nor_Doppler_refineSearch_set,N,L, ...
    L_delay,L_Doppler,max_norDelay,max_norDoppler,N_delay_peak,N_Doppler_peak,fs_up,delta_f, ...
    Num_iter_single,Num_iter_cyclic,Th_noTarget,Max_num_target)
%%  版本信息
%%
%%2025/8/3 V1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2025  BUAA
%======================================================
% Designer：Ni Yuanhan
% Design Date：202500803
% Modifier：Ni Yuanhan
% Modifying Date：202500803
% Function
% Description：基于梯度下降估计目标时延、多普勒和增益
% Programming
% Idea：
% Input parameter:
% Subordination：无
% Others：
%======================================================
% Modification：f_paraEstibyGradDesc_newtonOMP：适用单目标，用峰值估计增益
%               f_paraEstibyGradDesc_newtonOMP_v1：适用单目标，用矩阵范数估计增益
%               f_paraEstibyGradDesc_newtonOMP_multiTarget_v1：适用多目标，用矩阵范数估计增益
% Record：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Ts = 1/(fs_up);

    % 参数估计结果集合 (时延；多普勒；增益)
    set_nor_parameters_esti = [];

    % 记录每个目标峰值对应的二维索引
    set_index_delay_set_fracEst = [];
    set_index_Doppler_set_fracEst = [];

    [Len_x,Len_y] = size(ML_value_refine_set_aver);

    max_ML_value = max(max(ML_value_refine_set_aver));
    Th_peak_dB = 21;
    % 考虑到强目标消除后，后有残留误差，如果门限过低，会出现很多假目标
    k0 = 0.0044; % 经验值
    if Th_noTarget < k0*max_ML_value
        Th_noTarget = k0*max_ML_value;
    end

    % 启动迭代过程
    ML_value_residual = ML_value_refine_set_aver;
    num_outer = 1;
    while true
        num_outer = num_outer + 1;
        % 搜索峰值
        [value_row_max_refine, I_row_set_max_refine] = max(ML_value_residual);
        [value_max_refine, I_col_max_refine] = max(value_row_max_refine);
        I_row_max_refine = I_row_set_max_refine(I_col_max_refine);

        % CFAR 检测是否存在目标
        if value_max_refine <= Th_noTarget % 不存在目标
            break;
        end

%         figure
%         mesh(10*log10(abs(ML_value_residual)));

        % 计算时延和多普勒估计值
        nor_delay_est = Nor_delay_refineSearch_set(I_row_max_refine);
        nor_Doppler_est = Nor_Doppler_refineSearch_set(I_col_max_refine);
        % 增益粗估计值
        P_sig_strong_estiRough = value_max_refine;

        % 只用峰值左右 N_delay_peak 和 N_Doppler_peak 个点估计非整数时延和多普勒
    %         N_delay_peak = 2;
    %         N_Doppler_peak = 4;

        % 做一便异常值检测，防止索引越界(保证目标不处于异常区域内)

        if I_row_max_refine + N_delay_peak > Len_x
            I_row_max_refine = Len_x - N_delay_peak;
        elseif I_row_max_refine - N_delay_peak <= 0
            I_row_max_refine = N_delay_peak + 1;
        end
        if I_col_max_refine + N_Doppler_peak > Len_y
            I_col_max_refine = Len_y - N_Doppler_peak;
        elseif I_col_max_refine - N_Doppler_peak <= 0
            I_col_max_refine = N_Doppler_peak + 1;
        end
        % 索引
        index_delay_set_fracEst = I_row_max_refine-N_delay_peak :1: I_row_max_refine+N_delay_peak;
        index_Doppler_set_fracEst = I_col_max_refine-N_Doppler_peak :1: I_col_max_refine+N_Doppler_peak;
        %% ---------------------------------------------------

        %% *********** Single Refinement 利用梯度下降法估计 非整数时延和多普勒 ******************
        L_delay_noDelayDop = 2;% 离线生成模糊函数时 对应的 时延维 过采样倍数
        L_Doppler_noDelayDop = 2;% 离线生成模糊函数时 对应的 多普勒维 过采样倍数
        L_delay_aver = 1/(Nor_delay_refineSearch_set(2)-Nor_delay_refineSearch_set(1)); % 仿真对应的 时延维 过采样倍数
        L_Doppler_aver = 1/(Nor_Doppler_refineSearch_set(2)-Nor_Doppler_refineSearch_set(1)); % 仿真对应的 多普勒维 过采样倍数

        % 对应的归一化时延和多普勒值
        nor_tau_list = Nor_delay_refineSearch_set(index_delay_set_fracEst)/L_delay_aver*L_delay_noDelayDop;
        nor_v_list = Nor_Doppler_refineSearch_set(index_Doppler_set_fracEst)/L_Doppler_aver*L_Doppler_noDelayDop;

        % 对应上采样后的索引
    %     index_delay_set_fracEst_upSam = index_delay_set_fracEst/L_delay_aver*L_delay_noDelayDop;
    %     index_Doppler_set_fracEst_upSam = index_Doppler_set_fracEst/L_delay_aver*L_delay_noDelayDop;

        % 时延和多普勒间隔
        Delta_tau = 1/(fs_up)/L_delay_aver /(L_delay_noDelayDop/L_delay_aver);
        Delta_v = delta_f /L_Doppler_aver / (L_Doppler_noDelayDop/L_Doppler_aver);

        % 梯度函数
    %         Num_iter_single = Num_iter_set(index_PSWR);

        % 数据加噪
        D_func = ML_value_residual(index_delay_set_fracEst,index_Doppler_set_fracEst) ;
        % 参数初始化
        P_sig_strong_estiRefine_set = zeros(Num_iter_single,1);
        tau0_cal = zeros(Num_iter_single,1);
        v0_cal = zeros(Num_iter_single,1);
        tau0_cal(1) = nor_delay_est;
        v0_cal(1) = nor_Doppler_est;
        P_sig_strong_Refine_ini = 1;
        P_sig_strong_estiRefine_set(1) = P_sig_strong_estiRough;
        del_iter = zeros(Num_iter_single,1);
        num = 1;
        while true
            num = num + 1;

            if num > 25
               stop = 1;
            end

            [Ja_matrix, He_matrix] = f_grad2DAmbigFunc_dealyDoppler_noParfor_newton_v2(tau0_cal(num-1), v0_cal(num-1), AmbigFunc_pulsedAFDM_noDelayDop, D_func, nor_tau_list, nor_v_list,index_delay_set_fracEst,index_Doppler_set_fracEst,Ts,delta_f, Delta_tau, Delta_v,N,L,L_delay,L_Doppler,L_delay_noDelayDop,L_Doppler_noDelayDop,max_norDelay,max_norDoppler,P_sig_strong_estiRefine_set(num-1));

            inv_He_matrix = He_matrix^(-1);
            del_delayDoppler = inv_He_matrix * Ja_matrix;

            tau0_cal(num) = tau0_cal(num-1) - del_delayDoppler(1);
            v0_cal(num) = v0_cal(num-1) - del_delayDoppler(2);

            % 估计目标增益
            nor_delay_cal = tau0_cal(num);
            nor_Doppler_cal = v0_cal(num);
            % 对应过采样的结果
            nor_delay_esti_upSam = nor_delay_cal * L_delay_aver * (L_delay_noDelayDop/L_delay_aver);
            nor_Doppler_esti_upSam = nor_Doppler_cal * L_Doppler_aver * (L_Doppler_noDelayDop/L_Doppler_aver);
            % 反演强目标回波
            D_func_search_tmp = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor(AmbigFunc_pulsedAFDM_noDelayDop,N,L, ...
                L_delay_noDelayDop,L_Doppler_noDelayDop,P_sig_strong_Refine_ini,nor_delay_esti_upSam,nor_Doppler_esti_upSam);
            echo_strong_esti_tmp  = D_func_search_tmp(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);

            % 截取目标附近的区域，用于估计幅度（功率量纲）
            echo_strong_esti = echo_strong_esti_tmp(index_delay_set_fracEst,index_Doppler_set_fracEst);

            % figure;
            % [X, Y] = meshgrid(index_Doppler_set_fracEst, index_delay_set_fracEst);
            % % mesh(X, Y, D_func_migration);
            % mesh(X, Y, 10 * log10(echo_strong_esti / max(max(echo_strong_esti))));
            % ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
            % xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
            % zlabel('Normalized average squared MF');
            % set(gca, 'FontSize', 10);
            % % title('接收回波匹配滤波结果');

            P_sig_strong_estiRefine_set(num) = sum(sum(D_func .* echo_strong_esti))./norm(echo_strong_esti,'fro')^2;

            del_iter(num) = norm([tau0_cal(num),v0_cal(num),P_sig_strong_estiRefine_set(num)] - [tau0_cal(num-1),v0_cal(num-1),P_sig_strong_estiRefine_set(num-1)]);
            if del_iter(num) < 1e-10 || (num >= Num_iter_single)
                tau0_esti = tau0_cal(num);
                v0_esti = v0_cal(num);
                nor_delay_esti = tau0_esti;
                nor_Doppler_esti = v0_esti;
                P_sig_strong_estiRefine = P_sig_strong_estiRefine_set(num);
                amp_sig_strong_esti = sqrt(P_sig_strong_estiRefine);

                % 跟新参数估计集合
                set_nor_parameters_esti = [set_nor_parameters_esti, [nor_delay_esti;nor_Doppler_esti;P_sig_strong_estiRefine]];

                % 记录每个目标峰值对应的二维索引
                set_index_delay_set_fracEst = [set_index_delay_set_fracEst; index_delay_set_fracEst];
                set_index_Doppler_set_fracEst = [set_index_Doppler_set_fracEst; index_Doppler_set_fracEst];
                break;
            end
        end
        %% ================ Single Refinement ================================%

        %% **************** Cyclic Refinement ********************************
        Num_detected_target = size(set_nor_parameters_esti,2);
        num_cyclic = 0;
        while true
            num_cyclic = num_cyclic + 1;
            for n = Num_detected_target:-1:1
                % 第 n 个目标的参数
                nor_delay_esti_n = set_nor_parameters_esti(1,n);
                nor_Doppler_esti_n = set_nor_parameters_esti(2,n);
                P_esti_n = set_nor_parameters_esti(3,n);
                index_delay_set_fracEst_n = set_index_delay_set_fracEst(n,:);
                index_Doppler_set_fracEst_n = set_index_Doppler_set_fracEst(n,:);
                nor_tau_list_n = Nor_delay_refineSearch_set(index_delay_set_fracEst_n)/L_delay_aver*L_delay_noDelayDop;
                nor_v_list_n = Nor_Doppler_refineSearch_set(index_Doppler_set_fracEst_n)/L_Doppler_aver*L_Doppler_noDelayDop;

                % 反演出其他目标的匹配滤波结果
                echo_otherTarget_esti_set = zeros(Len_x,Len_y,Num_detected_target-1);
                k_num = 0;
                for k = 1:Num_detected_target
                    if k ~= n % 非第 n 个目标的参数
                        k_num = k_num + 1;
                        nor_delay_esti_k = set_nor_parameters_esti(1,k);
                        nor_Doppler_esti_k = set_nor_parameters_esti(2,k);
                        P_esti_k = set_nor_parameters_esti(3,k);
                        % 对应过采样的结果
                        nor_delay_esti_upSam = nor_delay_esti_k * L_delay_aver * (L_delay_noDelayDop/L_delay_aver);
                        nor_Doppler_esti_upSam = nor_Doppler_esti_k * L_Doppler_aver * (L_Doppler_noDelayDop/L_Doppler_aver);
                        % 反演强目标回波
                        D_func_search_tmp = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor(AmbigFunc_pulsedAFDM_noDelayDop,N,L,L_delay_noDelayDop,L_Doppler_noDelayDop,P_esti_k,nor_delay_esti_upSam,nor_Doppler_esti_upSam);
                        echo_strong_esti_tmp  = D_func_search_tmp(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
                        echo_otherTarget_esti_set(:,:,k_num) = echo_strong_esti_tmp;
                    end
                end

                % 用接收的匹配滤波结果，消除其他回波影响，得到第n个目标对应的接收的匹配滤波结果
                echo_Target_n_esti = (ML_value_refine_set_aver - sum(echo_otherTarget_esti_set,3));

%                 figure
%                 mesh(10*log10(abs(echo_Target_n_esti)));

                % 再次 Refine 第n个目标的参数
                % 数据加噪
                D_func = echo_Target_n_esti(index_delay_set_fracEst_n,index_Doppler_set_fracEst_n) ;
                % 参数初始化
                P_sig_strong_estiRefine_set = zeros(Num_iter_single,1);
                tau0_cal = zeros(Num_iter_single,1);
                v0_cal = zeros(Num_iter_single,1);
                tau0_cal(1) = nor_delay_esti_n;
                v0_cal(1) = nor_Doppler_esti_n;
                P_sig_strong_Refine_ini = 1;
                P_sig_strong_estiRefine_set(1) = P_esti_n;
                del_iter = zeros(Num_iter_single,1);
                num = 1;
                while true
                    num = num + 1;

                    if num > 25
                       stop = 1;
                    end

                    [Ja_matrix, He_matrix] = f_grad2DAmbigFunc_dealyDoppler_noParfor_newton_v2(tau0_cal(num-1), v0_cal(num-1), AmbigFunc_pulsedAFDM_noDelayDop, D_func, nor_tau_list_n, nor_v_list_n,index_delay_set_fracEst_n,index_Doppler_set_fracEst_n,Ts,delta_f, Delta_tau, Delta_v,N,L,L_delay,L_Doppler,L_delay_noDelayDop,L_Doppler_noDelayDop,max_norDelay,max_norDoppler,P_sig_strong_estiRefine_set(num-1));

                    inv_He_matrix = He_matrix^(-1);
                    del_delayDoppler = inv_He_matrix * Ja_matrix;

                    tau0_cal(num) = tau0_cal(num-1) - del_delayDoppler(1);
                    v0_cal(num) = v0_cal(num-1) - del_delayDoppler(2);

                    % 估计目标增益
                    nor_delay_cal = tau0_cal(num);
                    nor_Doppler_cal = v0_cal(num);
                    % 对应过采样的结果
                    nor_delay_esti_upSam = nor_delay_cal * L_delay_aver * (L_delay_noDelayDop/L_delay_aver);
                    nor_Doppler_esti_upSam = nor_Doppler_cal * L_Doppler_aver * (L_Doppler_noDelayDop/L_Doppler_aver);
                    % 反演强目标回波
                    D_func_search_tmp = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor(AmbigFunc_pulsedAFDM_noDelayDop,N,L,L_delay_noDelayDop,L_Doppler_noDelayDop,P_sig_strong_Refine_ini,nor_delay_esti_upSam,nor_Doppler_esti_upSam);
                    echo_strong_esti_tmp  = D_func_search_tmp(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);

                    % 截取目标附近的区域，用于估计幅度（功率量纲）
                    echo_strong_esti = echo_strong_esti_tmp(index_delay_set_fracEst_n,index_Doppler_set_fracEst_n);
                    P_sig_strong_estiRefine_set(num) = sum(sum(D_func .* echo_strong_esti))./norm(echo_strong_esti,'fro')^2;

                    del_iter(num) = norm([tau0_cal(num),v0_cal(num),P_sig_strong_estiRefine_set(num)] - [tau0_cal(num-1),v0_cal(num-1),P_sig_strong_estiRefine_set(num-1)]);
                    if del_iter(num) < 1e-10 || (num >= Num_iter_single)
%                     if del_iter(num) < 1e-10 || (num >= Num_iter_cyclic)
                        tau0_esti = tau0_cal(num);
                        v0_esti = v0_cal(num);
                        nor_delay_esti = tau0_esti;
                        nor_Doppler_esti = v0_esti;
                        P_sig_strong_estiRefine = P_sig_strong_estiRefine_set(num);
                        amp_sig_strong_esti = sqrt(P_sig_strong_estiRefine);

                        % 跟新第n个目标的参数
                        set_nor_parameters_esti(:,n) = [nor_delay_esti;nor_Doppler_esti;P_sig_strong_estiRefine];
                        break;
                    end
                end
            end
            if num_cyclic >= Num_iter_cyclic
               break;
            end
        end
        %% ================ Cyclic Refinement ================================%

        %% *************** SIC 所有检测到的目标影响，更新残余的匹配滤波结果 *********************%
        % 反演出所有被检测目标的匹配滤波结果
        echo_detectTarget_esti_set = zeros(Len_x,Len_y,Num_detected_target);
        for k = 1:Num_detected_target
            nor_delay_esti_k = set_nor_parameters_esti(1,k);
            nor_Doppler_esti_k = set_nor_parameters_esti(2,k);
            P_esti_k = set_nor_parameters_esti(3,k);
            % 对应过采样的结果
            nor_delay_esti_upSam = nor_delay_esti_k * L_delay_aver * (L_delay_noDelayDop/L_delay_aver);
            nor_Doppler_esti_upSam = nor_Doppler_esti_k * L_Doppler_aver * (L_Doppler_noDelayDop/L_Doppler_aver);
            % 反演强目标回波
            D_func_search_tmp = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor(AmbigFunc_pulsedAFDM_noDelayDop,N,L,L_delay_noDelayDop,L_Doppler_noDelayDop,P_esti_k,nor_delay_esti_upSam,nor_Doppler_esti_upSam);
            echo_strong_esti_tmp  = D_func_search_tmp(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
            echo_detectTarget_esti_set(:,:,k) = echo_strong_esti_tmp;
        end
        % 用接收的匹配滤波结果，消除所有被检测目标，更新残余的匹配滤波结果
        ML_value_residual = (ML_value_refine_set_aver - sum(echo_detectTarget_esti_set,3));

%         figure
%         mesh((abs(ML_value_residual)));

%         del_outer = norm(ML_value_residual,'fro')^2/(Len_x*Len_y);
        max_ML_value_residual = max(max(ML_value_residual));
        if max_ML_value_residual <= Th_noTarget || Num_detected_target >= Max_num_target
            break;
        end
    end
%             % 估计目标增益
%             % 对应过采样的结果
%             nor_delay_esti_upSam = nor_delay_esti * L_delay_aver * (L_delay_noDelayDop/L_delay_aver);
%             nor_Doppler_esti_upSam = nor_Doppler_esti * L_Doppler_aver * (L_Doppler_noDelayDop/L_Doppler_aver);
%             % 反演强目标回波
%             D_func_search_tmp = f_grad2DAmbigFunc_fracDelayFracDoppler_value(AmbigFunc_pulsedAFDM_noDelayDop,N,L,L_delay_noDelayDop,L_Doppler_noDelayDop,P_sig_strong_Refine_ini,nor_delay_esti_upSam,nor_Doppler_esti_upSam);
%             echo_strong_esti_tmp  = D_func_search_tmp(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
%             % 读取峰值处的幅度（功率量纲）
%             amp_echo_strong = echo_strong_esti_tmp(I_row_max_refine,I_col_max_refine);
%             P_sig_strong_estiRefine_set(num_amp) = P_sig_strong_Refine_ini/amp_echo_strong;
%
%             del_iter_amp = abs(P_sig_strong_estiRefine_set(num_amp) - P_sig_strong_estiRefine_set(num_amp-1))/P_sig_strong_Refine_ini;
%             if del_iter_amp < 1e-6 || (num_amp >= Num_iter_amp)
%                 P_sig_strong_estiRefine = P_sig_strong_estiRefine_set(num_amp);
%                 break;
%             end
%         end
        % 目标增益的估计值
%         P_sig_strong_esti = P_sig_strong_estiRough * P_sig_strong_estiRefine;
%         amp_sig_strong_esti = sqrt(P_sig_strong_estiRough * P_sig_strong_estiRefine);


%         figure
%         mesh(echo_strong_esti_tmp);
%         figure
%         mesh(ML_value_refine_set_aver);
%         P_sig_strong_estiRefine_2 = sum(sum(ML_value_refine_set_aver .* echo_strong_esti_tmp))./norm(echo_strong_esti_tmp,'fro')^2;
%         amp_sig_strong_esti_2 = sqrt(P_sig_strong_estiRefine_2);
%                 sum(sum(echo_strong_esti_tmp .* echo_strong_esti_tmp))./norm(echo_strong_esti_tmp,'fro')^2;


end