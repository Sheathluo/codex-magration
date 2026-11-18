%%% 2026年08月23日 创建人：罗韬
%%%% 根据时延和多普勒画出损失值

% true_delay = 5.22346;
% delay_set = true_delay + (-1:0.1:1);
% true_dopple = 5.85734;
% doppler_set = true_dopple + (-1:0.1:1);
%
% delay_len = length(delay_set);
% doppler_len = length(doppler_set);
% cost_vec = zeros([delay_len,doppler_len]);
% P_vec = zeros([delay_len,doppler_len]);
%
% parfor delay_index = 1:delay_len
%     for dopple_index = 1:doppler_len
%         tau_new_migration = delay_set(delay_index) + (0:para.Nsym-1).' .* gamma * doppler_set(dopple_index);
%
%         nor_delay_esti_upSam = tau_new_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
%         nor_Doppler_esti_upSam = doppler_set(dopple_index) * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);
%
%         echo_unit = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v2(AmbigFunc_pulsedAFDM_noDelayDop, ...
%             1, nor_delay_esti_upSam, nor_Doppler_esti_upSam, para, echo_rows, echo_cols, cfg);
%
%         echo_local = echo_unit(index_delay_set_fracEst, index_Doppler_set_fracEst);
%         P_new = sum(D_func .* echo_local, "all") ./ norm(echo_local, 'fro')^2;
%         P_vec(dopple_index, delay_index) = P_new;
%         cost = norm(D_func - P_new * echo_local, 'fro')^2;
%         cost_vec(delay_index, dopple_index) = cost;
%     end
% end
%
% % figure
% % plot(delay_set, cost_vec, 'LineWidth', 1.6);
% % hold on;
% % xline(true_delay, '--r', '真实时延', 'LineWidth', 1.2);
% %
% % plot(doppler_set, cost_vec, 'LineWidth', 1.6);
% % hold on;
% % xline(true_dopple, '--b', '真实多普勒', 'LineWidth', 1.2);
%
% figure
% mesh(delay_set, doppler_set, cost_vec)



true_delay   = 5.22346;
delay_set    = true_delay + (-1:0.1:1);
true_dopple  = 5.85734;
doppler_set  = true_dopple + (-1:0.1:1);

delay_len   = length(delay_set);
doppler_len = length(doppler_set);
cost_vec = zeros(delay_len, doppler_len);
P_vec    = zeros(delay_len, doppler_len);

parfor delay_index = 1:delay_len
    for dopple_index = 1:doppler_len
        tau_new_migration = delay_set(delay_index) + (0:para.Nsym-1).' .* gamma * doppler_set(dopple_index);

        nor_delay_esti_upSam   = tau_new_migration * L_delay_aver * (L_delay_noDelayDop / L_delay_aver);
        nor_Doppler_esti_upSam = doppler_set(dopple_index) * L_Doppler_aver * (L_Doppler_noDelayDop / L_Doppler_aver);

        echo_unit = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_v2( ...
            AmbigFunc_pulsedAFDM_noDelayDop, 1, nor_delay_esti_upSam, ...
            nor_Doppler_esti_upSam, para, echo_rows, echo_cols, cfg);

        echo_local = echo_unit(index_delay_set_fracEst, index_Doppler_set_fracEst);
        P_new = sum(D_func .* echo_local, "all") ./ norm(echo_local, 'fro')^2;

        P_vec(delay_index, dopple_index)    = P_new;                       % 修正索引
        cost_vec(delay_index, dopple_index) = norm(D_func - P_new * echo_local, 'fro')^2;
    end
end

%% 三维代价函数图
figure;
mesh(delay_set, doppler_set, cost_vec.');   % 转置，使 Z 为 doppler_len×delay_len
hold on;

% 真实值对应的代价（取网格上离真实值最近的点）
[~, idx_d] = min(abs(delay_set   - true_delay));
[~, idx_p] = min(abs(doppler_set - true_dopple));
cost_true  = cost_vec(idx_d, idx_p);

% 在三维图上标出真实值
plot3(true_delay, true_dopple, cost_true, 'rp', ...
    'MarkerSize', 14, 'MarkerFaceColor', 'r');
text(true_delay, true_dopple, cost_true * 1.05, ...
    sprintf('真值 (%.4f, %.4f)', true_delay, true_dopple), ...
    'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');

xlabel('时延 (Delay)', 'FontSize', 12);
ylabel('多普勒 (Doppler)', 'FontSize', 12);
zlabel('代价函数 Cost', 'FontSize', 12);
title('二维时延-多普勒代价函数', 'FontSize', 13);
grid on;
view(45, 30);
colorbar;