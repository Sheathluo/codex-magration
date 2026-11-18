function E_mat = f_theoretical_AmbigFunc_pulseShaped_AFDM_v6(p, N, L, c1, mu4, tau_list, v_list, tau0, v0)
% v5:在v3的基础上，将线性移位改成循环移位，将对接收脉冲移位改成对参考脉冲移位
% v6:在v5的基础上，将tau0和v0作用在接收信号上，具备生成回波匹配滤波理论结果的功能
% compute_EX2 计算公式 E{|X(tau,v)|^2}
% 输入:
%   p    - 长度为 NL 的实数脉冲成形向量
%   N    - 符号数
%   L    - 参数 L
%   c1   - 参数 c1
%   mu4  - 参数 mu4
%   tau  - 时延 tau
%   v    - 多普勒 v
% 输出:
%   E_X_tau_v_sq - 结果

    NL = N * L;
    n_vec = 0:N-1;
    m_vec = 0:N-1;

%     tau_list = 0;
%     v_list = 0;

    % 非整数时延 sinc(n)=sin[pi (n-x)]/[pi (n- x)]函数的抽头数量 2*L_filter + 1
    L_filter = 5000;
%         L_filter = 50000;
    coe_sinc = zeros(2*L_filter+1,1);
    n_set = -L_filter:1:L_filter;


%     figure
%     plot(real(p),'-b*'); hold on;
%     plot(real(p_Tx),'--ro'); hold on;
%     grid on;


    % 初始化结果矩阵
    E_mat_q = zeros(length(tau_list), length(v_list));

    for ti = 1:length(tau_list)
        ti
        tau = tau_list(ti);

%         int_tau = round(tau);
%         frac_tau = tau - int_tau;

%         if abs(frac_tau) < 1e-10
%             coe_sinc(:,1) = [zeros(L_filter,1);1;zeros(L_filter,1)];
%         else
%             for i = 1:2*L_filter+1
%                 n = n_set(i);
% %                 coe_sinc(i,1) = sin(pi*(n-frac_tau))/(pi*(n-frac_tau));
%                 coe_sinc(i,1) = sin(pi*(n+frac_tau))/(pi*(n+frac_tau));
%     %                 coe_sinc(i,1) = abs(sin(pi*(n-frac_norDelay))/(pi*(n-frac_norDelay)));
%     %                 coe_sinc(i,1) = sin(pi*(n-frac_norDelay))/sin(pi/(N*L)*(n-frac_norDelay));
%             end
%         end
%         coe_sinc = coe_sinc/norm(coe_sinc);
%
%         p_Tx_tmp = conv(coe_sinc,fftshift(p));
%         p_Tx = fftshift(p_Tx_tmp(L_filter+1:L_filter+length(p)));

        % p_Tx：表示发送信号到接收端的回波信号，直接用fft形成循环移位
        p_Tx = real( ifft(fftshift(fftshift(fft(p)).*exp(-1i*2*pi./(length(p)).*(tau0)*((0:length(p)-1)-length(p)/2).'))));

        % p_Rx：表示接收端匹配滤波的参考信号，直接用fft形成循环移位
        p_Rx = real( ifft(fftshift(fftshift(fft(p)).*exp(-1i*2*pi./(length(p)).*(tau)*((0:length(p)-1)-length(p)/2).'))));


%         % 线性移位和循环移位的差异
%         del_p_Tx_12 = norm(p_Tx - p_Tx_2);

%         figure
%         plot(real(p),'-r*'); hold on;
%         plot(real(p_Tx),'--bo'); hold on;
%         plot(real(p_Tx_2),'-.ks'); hold on;
%         grid on;
%         legend('原响应','线性移位后响应','循环移位后响应');


%         for vi = 1:length(v_list)
        parfor vi = 1:length(v_list)
            v = v_list(vi);
            p_Tx_tmp = p_Tx;
            p_Rx_tmp = p_Rx;

            % ------------------------
            % Part 1
            sum1 = 0;
            for n = n_vec
                inner_sum = 0;
                for k = -n*L : NL - n*L - 1
                    k_mod = mod(k, NL) + 1;
                    ktau_mod = mod(k - 0, NL) + 1;
%                     inner_sum = inner_sum + conj(p_Tx(k_mod)) * p(ktau_mod) * exp(1j*2*pi/NL*v*k); % 对接收脉冲移位
%                     inner_sum = inner_sum + conj(p_tmp(k_mod)) * p_Rx_tmp(ktau_mod) * exp(1j*2*pi/NL*v*k); % 对参考脉冲移位
                    inner_sum = inner_sum + (p_Tx_tmp(k_mod)) * conj(p_Rx_tmp(ktau_mod)) * exp(-1j*2*pi/NL*(v-v0)*k); % tau_0 对接收脉冲移位, tau对参考脉冲移位
                end
                sum1 = sum1 + exp(-1j * 2*pi/N * (v-v0) * n) * inner_sum;
            end
            part1 = abs(sum1)^2;

            % ------------------------
            % Part 2
%             part2_tmp = 0;
%             for q = m_vec
%                 for q1 = m_vec
%                     sum_nm = 0;
%                     for n = n_vec
%                         exponent = exp(-1j*2*pi/N*(2*N*c1*(n-q) - v)*n);
%                         inner_sum = 0;
%                         for k = -n*L : NL - n*L - 1
%                             k_mod = mod(k, NL) + 1;
%                             ktau_mod = mod(k - (tau - (n-q)*L), NL) + 1;
%                             inner_sum = inner_sum + conj(p(k_mod)) * p(ktau_mod) * exp(1j*2*pi/NL*v*k);
%                         end
%                         sum_nm = sum_nm + exponent * inner_sum;
%                     end
%                     part2_tmp = part2_tmp + abs(sum_nm)^2;
%                 end
%             end
            part2_tmp = 0;
            for m = m_vec
                sum_nm = 0;
                for n = n_vec
                    exponent = exp(1j*2*pi/N*(2*N*c1*(m) - (v-v0))*n);
                    inner_sum = 0;
                    for k = -n*L : NL - n*L - 1
                        k_mod = mod(k, NL) + 1;
                        ktau_mod = mod(k - (0 - (m)*L), NL) + 1;
%                         inner_sum = inner_sum + conj(p_Tx(k_mod)) * p(ktau_mod) * exp(1j*2*pi/NL*v*k); % 对接收脉冲移位
%                         inner_sum = inner_sum + conj(p_tmp(k_mod)) * p_Rx_tmp(ktau_mod) * exp(1j*2*pi/NL*v*k); % 对参考脉冲移位
                        inner_sum = inner_sum + (p_Tx_tmp(k_mod)) * conj(p_Rx_tmp(ktau_mod)) * exp(-1j*2*pi/NL*(v-v0)*k); % tau_0 对接收脉冲移位, tau对参考脉冲移位
                    end
                    sum_nm = sum_nm + exponent * inner_sum;
                end
                part2_tmp = part2_tmp + abs(sum_nm)^2;
            end
            part2 = part2_tmp/ N * (mu4 - 2) ;

            % ------------------------
            % Part 3
            part3 = 0;
            for m = m_vec
                for n = n_vec
                    inner_sum = 0;
                    for k = -n*L : NL - n*L - 1
                        k_mod = mod(k, NL) + 1;
                        ktau_mod = mod(k - (0 - (m)*L), NL) + 1;
%                         inner_sum = inner_sum + conj(p_Tx(k_mod)) * p(ktau_mod) * exp(1j*2*pi/NL*v*k); % 对接收脉冲移位
%                         inner_sum = inner_sum + conj(p_tmp(k_mod)) * p_Rx_tmp(ktau_mod) * exp(1j*2*pi/NL*v*k); % 对参考脉冲移位
                        inner_sum = inner_sum + (p_Tx_tmp(k_mod)) * conj(p_Rx_tmp(ktau_mod)) * exp(-1j*2*pi/NL*(v-v0)*k); % tau_0 对接收脉冲移位, tau对参考脉冲移位
                    end
                    part3 = part3 + abs(inner_sum)^2;
                end
            end

            % ------------------------
            % 存储结果
            E_mat_q(ti, vi) = part1 + part2 + part3;
        end
    end
    E_mat(:,:) = E_mat_q;

end
