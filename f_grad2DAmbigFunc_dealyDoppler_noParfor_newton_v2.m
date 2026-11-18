function [J, H] = f_grad2DAmbigFunc_dealyDoppler_noParfor_newton_v2(tau0, v0, X, Y, nor_tau_list, nor_v_list,...
    index_delay_set_fracEst,index_Doppelr_set_fracEst,Ts,delta_f, Delta_tau, ...
    Delta_v, N,L,L_delay,L_Doppler,L_delay_noDelayDop,L_Doppler_noDelayDop, ...
    max_norDelay,max_norDoppler,P_sig_strong_nor)
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
% Description：验证根据推导的模糊函数解析式，求梯度。适用非相干累积
% Programming
% Idea：
% Input parameter:
% Subordination：无
% Others：
%======================================================
% Modification：在 f_grad2DAmbigFunc_dealyDoppler_noParfor_newton 的基础上，
% 将迭代估计绝对时延和多普勒，改为估计归一化的时延和多普勒，看能否解决Hession矩阵求逆出奇异值的问题
% Record：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 计算 f(tau0, v0) 以及梯度 ∂f/∂tau0, ∂f/∂v0
% X, Y 是 N_delay×N_Doppler 矩阵
% tau: N×1 向量，v: 1×M 向量
% Delta_tau, Delta_v: 步长

[N_delay, N_Doppler] = size(X);
[T,V] = size(Y);
% Z = zeros(T, V);
dZ_dtau0 = zeros(T, V);
dZ_dv0 = zeros(T, V);

epsilon = 1e-8; % 奇点检测
% 相对于在时延和多普勒域过采样后的归一化时延和多普勒
nor_tau0 = tau0/Delta_tau;
nor_v0 = v0/Delta_v;

delta_f = Delta_v * L_Doppler_noDelayDop;

% 相对时延和多普勒分辨率归一化的时延和多普勒
nor_Doppler = v0/Delta_v/L_Doppler_noDelayDop;

% D_func_search_tmp0 = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor(X,N,L,L_delay_noDelayDop,L_Doppler_noDelayDop,1,nor_tau0,nor_v0);
% % D_func_search_tmp0 = f_grad2DAmbigFunc_fracDelayFracDoppler_value_noParfor_coherent(X,N,Ncp,Nsym,N_preSuffix,L,L_delay_noDelayDop,L_Doppler_noDelayDop,P_sig_strong_nor,nor_tau0,nor_v0);
% % D_func_search_tmp  = D_func_search_tmp0(1:L_delay_noDelayDop:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
% D_func_search_tmp  = D_func_search_tmp0(1:L_delay_noDelayDop/L_delay:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
% Z_1_0 = D_func_search_tmp(index_delay_set_fracEst,index_Doppelr_set_fracEst);
%
% Z_0 = zeros(size(Y,1),size(Y,2),size(Y,3));
% Z_2_0 = zeros(1,size(Y,3));
% for index_page = 1:length(index_DoppelrDash_set_fracEst)
% %                     theta = theta_peak - del_DoppelrDash_set(index_page);
%     theta = Nsym/N*nor_Doppler*(N+Ncp+2*N_preSuffix) - (index_DoppelrDash_set_fracEst(index_page)-1) - Nsym/2;
%     if mod(theta, Nsym) == 0
%         Z_2_0(1,index_page) = 1;
%     else
%         Z_2_0(1,index_page) = abs(1/Nsym*(exp(1i*2*pi*theta)-1)/(exp(1i*2*pi/Nsym*theta)-1))^2;
%     end
%     Z_0(:,:,index_page) = Z_1_0 * Z_2_0(1,index_page);
% end

% 生成傅里叶矩阵
F_tau = 1/sqrt(N*L*L_delay_noDelayDop) * fft(eye(N*L*L_delay_noDelayDop));  % 生成 N×N 傅里叶矩阵
F_v = 1/sqrt(N*L*L_Doppler_noDelayDop) * fft(eye(N*L*L_Doppler_noDelayDop));  % 生成 N×N 傅里叶矩阵

%% ------------------- Z1 -------------------
% Lambda_tau
% tau0 = 0.25;
% m_vec_tau = [0:N*L*L_delay_noDelayDop-1];
m_vec_tau = [0:0.5*N*L*L_delay_noDelayDop-1, (-0.5*N*L*L_delay_noDelayDop):-1];
vec_Lambda_tau = exp(-1j*2*pi/(N*L*L_delay_noDelayDop) * tau0 * Ts/Delta_tau * m_vec_tau);
maxtrix_Lambda_tau_kron = kron(vec_Lambda_tau.',ones(1,size(X,2)));
Lambda_tau = diag(vec_Lambda_tau);


% Lambda_v
% m_vec_v = [0:N*L*L_Doppler_noDelayDop-1];
m_vec_v = [0:0.5*N*L*L_Doppler_noDelayDop-1, (-0.5*N*L*L_Doppler_noDelayDop):-1];
vec_Lambda_v = exp(-1j*2*pi/(N*L*L_Doppler_noDelayDop) * v0 * delta_f/Delta_v * m_vec_v);
maxtrix_Lambda_v_kron = kron(ones(size(X,1),1),vec_Lambda_v);
Lambda_v = diag(vec_Lambda_v);

% Z1
Z1 = real( ifft(fft(ifft(maxtrix_Lambda_tau_kron .* fft(X,[],1), [],1),[],2) .* maxtrix_Lambda_v_kron,[],2));
% Z1_tmp = real( F_tau' * Lambda_tau * F_tau * X * F_v * Lambda_v * F_v' );
% del_Z1 = norm(Z1 - Z1_tmp);

% 验证频域乘以线性相位与时域循环卷积 产生的非整数循环移位是否一致
% % 频域
% X_delay_freq = F_tau' * Lambda_tau * F_tau * X;
%
% % 时域
% N_total = N*L*L_Doppler_noDelayDop;
% h_cir = zeros(N_total,1);
% for i=0:N_total-1
%     theta = 2*pi*(i - tau0 * Ts/Delta_tau);
%     if mod(theta, N_total) == 0
%         h_cir(i+1,1) = 1;
%     else
%         h_cir(i+1,1) = 1/N_total*(exp(1i*theta/2) - exp(-1i*theta/2))/(exp(1i*theta/N_total) - 1);
%     end
% end
% X_delay_time = zeros(N_total,1);
% for k = 0:N_total-1
%     tmp = zeros(N_total,1);
%     for m = 0:N_total-1
%         tmp(m+1,1) = X(m+1,1) .* h_cir(mod(k-m,N_total)+1,1);
%     end
%     X_delay_time(k+1,1) = sum(tmp);
% end
%
% del_Z1_tmp = norm(X_delay_time - X_delay_freq(:,1));
%
% figure
% plot(real(X_delay_freq(:,1)),'-r*'); hold on;
% plot(real(X_delay_time(:,1)),'--bo'); hold on;
% legend('频域实现','时域实现');
% -----------------------------------------------------------------%

% 一阶导 Z1
% dLambda_tau = diag(-1j*2*pi/(N*L*L_delay_noDelayDop) * m_vec_tau / Delta_tau .* exp(-1j*2*pi/(N*L*L_delay_noDelayDop) * tau0/Delta_tau * m_vec_tau));
% dLambda_v   = diag(-1j*2*pi/(N*L*L_Doppler_noDelayDop) * m_vec_v / Delta_v .* exp(-1j*2*pi/(N*L*L_Doppler_noDelayDop) * v0/Delta_v * m_vec_v));
% dZ1_dtau_0 = real(F_tau' * dLambda_tau * F_tau * X * F_v * Lambda_v * F_v');
% dZ1_dv_0   = real(F_tau' * Lambda_tau * F_tau * X * F_v * dLambda_v * F_v');

vec_dLambda_tau = -1j*2*pi/(N*L*L_delay_noDelayDop)*Ts/Delta_tau * m_vec_tau  .* exp(-1j*2*pi/(N*L*L_delay_noDelayDop) * tau0*Ts/Delta_tau * m_vec_tau);
maxtrix_dLambda_tau_kron = kron(vec_dLambda_tau.',ones(1,size(X,2)));
vec_dLambda_v = -1j*2*pi/(N*L*L_Doppler_noDelayDop)* delta_f/Delta_v * m_vec_v .* exp(-1j*2*pi/(N*L*L_Doppler_noDelayDop) * v0* delta_f/Delta_v * m_vec_v);
maxtrix_dLambda_v_kron = kron(ones(size(X,1),1),vec_dLambda_v);
dZ1_dtau = real( ifft(fft(ifft(maxtrix_dLambda_tau_kron .* fft(X,[],1), [],1),[],2) .* maxtrix_Lambda_v_kron,[],2));
dZ1_dv = real( ifft(fft(ifft(maxtrix_Lambda_tau_kron .* fft(X,[],1), [],1),[],2) .* maxtrix_dLambda_v_kron,[],2));


% 二阶导 Z1
% d2Lambda_tau = diag(-(2*pi/(N*L*L_delay_noDelayDop))^2 * (m_vec_tau/Delta_tau).^2 .* exp(-1j*2*pi/(N*L*L_delay_noDelayDop) * tau0/Delta_tau * m_vec_tau));
% d2Lambda_v   = diag(-(2*pi/(N*L*L_delay_noDelayDop))^2 * (m_vec_v/Delta_v).^2 .* exp(-1j*2*pi/(N*L*L_delay_noDelayDop) * v0/Delta_v * m_vec_v));
% d2Z1_dtau2_0 = real(F_tau' * d2Lambda_tau * F_tau * X * F_v * Lambda_v * F_v');
% d2Z1_dv2_0   = real(F_tau' * Lambda_tau * F_tau * X * F_v * d2Lambda_v * F_v');
% d2Z1_dtau_dv_0 = real(F_tau' * dLambda_tau * F_tau * X * F_v * dLambda_v * F_v');

vec_d2Lambda_tau = -(2*pi/(N*L*L_delay_noDelayDop))^2 * (m_vec_tau*Ts/Delta_tau).^2 .* exp(-1j*2*pi/(N*L*L_delay_noDelayDop) * tau0*Ts/Delta_tau * m_vec_tau);
maxtrix_d2Lambda_tau_kron = kron(vec_d2Lambda_tau.',ones(1,size(X,2)));
vec_d2Lambda_v = -1*(2*pi/(N*L*L_delay_noDelayDop))^2 * (m_vec_v* delta_f/Delta_v).^2 .* exp(-1j*2*pi/(N*L*L_delay_noDelayDop) * v0* delta_f/Delta_v * m_vec_v);
maxtrix_d2Lambda_v_kron = kron(ones(size(X,1),1),vec_d2Lambda_v);

d2Z1_dtau2 = real( ifft(fft(ifft(maxtrix_d2Lambda_tau_kron .* fft(X,[],1), [],1),[],2) .* maxtrix_Lambda_v_kron,[],2));
d2Z1_dv2 = real( ifft(fft(ifft(maxtrix_Lambda_tau_kron .* fft(X,[],1), [],1),[],2) .* maxtrix_d2Lambda_v_kron,[],2));
d2Z1_dtau_dv = real( ifft(fft(ifft(maxtrix_dLambda_tau_kron .* fft(X,[],1), [],1),[],2) .* maxtrix_dLambda_v_kron,[],2));


% 截取有效区间的数据
Z1  = Z1(1:L_delay_noDelayDop/L_delay:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
Z1 = Z1(index_delay_set_fracEst,index_Doppelr_set_fracEst);

dZ1_dtau  = dZ1_dtau(1:L_delay_noDelayDop/L_delay:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
dZ1_dtau = dZ1_dtau(index_delay_set_fracEst,index_Doppelr_set_fracEst);

dZ1_dv  = dZ1_dv(1:L_delay_noDelayDop/L_delay:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
dZ1_dv = dZ1_dv(index_delay_set_fracEst,index_Doppelr_set_fracEst);

d2Z1_dtau2  = d2Z1_dtau2(1:L_delay_noDelayDop/L_delay:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
d2Z1_dtau2 = d2Z1_dtau2(index_delay_set_fracEst,index_Doppelr_set_fracEst);

d2Z1_dv2  = d2Z1_dv2(1:L_delay_noDelayDop/L_delay:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
d2Z1_dv2 = d2Z1_dv2(index_delay_set_fracEst,index_Doppelr_set_fracEst);

d2Z1_dtau_dv  = d2Z1_dtau_dv(1:L_delay_noDelayDop/L_delay:max_norDelay*L*L_delay_noDelayDop,N*L*L_Doppler_noDelayDop/2-max_norDoppler*L*L_Doppler_noDelayDop/2+1:N*L*L_Doppler_noDelayDop/2+max_norDoppler*L*L_Doppler_noDelayDop/2);
d2Z1_dtau_dv = d2Z1_dtau_dv(index_delay_set_fracEst,index_Doppelr_set_fracEst);



%% ------------------- Z 及 一阶导 和 二阶导 -------------------
% Z = zeros(size(Y,1),size(Y,2));
% dZ_dtau = zeros(size(Y,1),size(Y,2));
% dZ_dv = zeros(size(Y,1),size(Y,2));
% d2Z_dtau2 = zeros(size(Y,1),size(Y,2));
% d2Z_dv2 = zeros(size(Y,1),size(Y,2));
% d2Z_dtau_dv = zeros(size(Y,1),size(Y,2));
% for index_m = 1: size(Y,3)
%     Z(:,:,index_m) = Z1 * Z2(1,index_m);
%     dZ_dtau(:,:,index_m) = Z2(1,index_m) * dZ1_dtau;
%     dZ_dv(:,:,index_m)   = Z2(1,index_m) * dZ1_dv + Z1 * dZ2_dv(1,index_m);
%
%     d2Z_dtau2(:,:,index_m) = Z2(1,index_m) * d2Z1_dtau2;
%     d2Z_dv2(:,:,index_m) = Z2(1,index_m) * d2Z1_dv2 + 2*dZ1_dtau*dZ1_dv + Z1*d2Z2_dv2(1,index_m);
%     d2Z_dtau_dv(:,:,index_m) = Z2(1,index_m) * d2Z1_dtau_dv + d2Z2_dv2(1,index_m)*dZ1_dtau;
% end

%% ************* 计算梯度 grad_tau0 和 grad_v0 *************
diff = Y - P_sig_strong_nor * Z1;

% figure;
% [X, Y1] = meshgrid(index_Doppelr_set_fracEst, index_delay_set_fracEst);
% mesh(X, Y1, Y);
% ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
% xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
% zlabel('Normalized average squared MF');
% % axis([-inf inf -inf inf -40 2]);
% set(gca, 'FontSize', 10);
% title('接收回波匹配滤波结果');
%
% figure;
% [X, Y1] = meshgrid(index_Doppelr_set_fracEst, index_delay_set_fracEst);
% mesh(X, Y1, Z1);
% ylabel('Normalized delay $\tau$', 'Interpreter', 'latex');
% xlabel('Normalized Doppler $\nu$', 'Interpreter', 'latex');
% zlabel('Normalized average squared MF');
% % axis([-inf inf -inf inf -40 2]);
% set(gca, 'FontSize', 10);
% title('理论模糊函数');

fval = sum(sum(sum(abs(diff))));
dS_dtau = 2*sum(sum(sum(diff .* dZ1_dtau * P_sig_strong_nor)));
dS_dv   = 2*sum(sum(sum(diff .* dZ1_dv * P_sig_strong_nor)));

d2S_dtau2 = -2*(sum(sum(sum( P_sig_strong_nor^2 * dZ1_dtau .* dZ1_dtau ))) - sum(sum(sum( - diff .* d2Z1_dtau2 * P_sig_strong_nor ))));
d2S_dv2 = -2*(sum(sum(sum( P_sig_strong_nor^2 * dZ1_dv .* dZ1_dv ))) - sum(sum(sum( - diff .* d2Z1_dv2 * P_sig_strong_nor ))));
d2S_dtau_dv = -2*(sum(sum(sum( P_sig_strong_nor^2 * dZ1_dtau .* dZ1_dv))) - sum(sum(sum(diff .* d2Z1_dtau_dv * P_sig_strong_nor ))));

% 雅可比矩阵
J = [dS_dtau; dS_dv];
% 海森矩阵
H = [d2S_dtau2, d2S_dtau_dv; d2S_dtau_dv, d2S_dv2];

end
