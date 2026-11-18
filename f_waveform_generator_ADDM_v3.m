function [try_catch, transmit_data_noCP, transmit_data_addCP, Tx_data] = f_waveform_generator_ADDM_v3(waveform_index, c1, c2, M, N, Ncp,Nsym, transmit_data)
%%  版本信息
%%
%%2024/04/08 V2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2024  BUAA
%======================================================
% Designer：Ni Yuanhan
% Design Date：20240408
% Modifier：Ni Yuanhan
% Modifying Date：20240527
% Function
% Description：在ADDM的基础上，支持生成全CP OTFS波形
% Programming
% Idea：
% Input parameter:  waveform_index - 波形类型
%                   c1 - 波形参数, 
%                   c2 - 波形参数，
%                   M - 每列有效符号个数
%                   N - 子载波个数
%                   Ncp - CP点数
%                   Nsym - 列数
%                   transmit_data - 发送数据
% Output parameter: try_catch - 错误指示标识, 
%                   transmit_data_noCP - 未添加CP的符号快, 
%                   transmit_data_addCP - 添加CP的符号快, 
%                   Tx_data - 生成的符号向量
% Subordination：无
% Others：
%======================================================
% Modification： 在v2的基础上，增加了生成AFDM符号块的功能
% Record：
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    try_catch = 0;
    if strcmp(waveform_index, 'ADDM')
        
    elseif strcmp(waveform_index, 'CPOTFS')
            
    elseif strcmp(waveform_index, 'AFDM')
%         c1 = (2*(alpha_max+xi_v)+1)/(2*N);
%         c2 = 1/(1000000*2*N);
    elseif strcmp(waveform_index, 'OFDM')
%         c1 = 0;
%         c2 = 0;
    else
        output = '请选择波形'
        Tx_data = 0;
        Rx_data = 0;
        Theta_c1 = 0;
        Theta_c2 = 0;
        try_catch = 1;
        return; % 跳出函数
    end
    Theta_c1 = zeros(N,N);
    Theta_c2 = zeros(M,M);
    F = zeros(M,N);
    for n=1:N
        Theta_c1(n,n) = exp(-1i*2*pi*c1*(n-1)^2);       
    end
    for m=1:M    
        Theta_c2(m,m) = exp(-1i*2*pi*c2*(m-1)^2);    
    end

%     for m=1:M
%        for n=1:N
%            F(m,n) = 1/sqrt(M)*exp(-1i*2*pi*(m-1)*(n-1)/M);
%        end
%     end

    % DAFT 变换矩阵
    % tic
    % A_0 = Theta_c2*F*Theta_c1; % 直接乘复杂度太高
    % toc

    % aa = diag(Theta_c1).';
%     kron_c1 = kron(ones(N,1),diag(Theta_c1).');
    % cc = diag(Theta_c2);
%     kron_c2 = kron(ones(1,M),diag(Theta_c2));
    % A = kron_c2.*F.*kron_c1; % 低复杂度计算

%     P_1 = var(transmit_data);

    %% **************** 开始生成ADDM波形 ************************
    % 首先二维符号逐行IDFT
    if strcmp(waveform_index, 'CPOTFS') || strcmp(waveform_index, 'ADDM')
        transmit_data_F = sqrt(Nsym)*ifft(transmit_data, Nsym, 2);
    else
        transmit_data_F = transmit_data;
    end
    
    if strcmp(waveform_index, 'CPOTFS')
        transmit_data_noCP = transmit_data_F; % 全CP OTFS
    else % AFDM/OFDM/OCDM
        % 然后逐列IDAFT
        kron_c2 = kron(ones(1,Nsym),diag(Theta_c2'));
        transmit_data_c2 = kron_c2.*transmit_data_F;% 低复杂度计算
        ifft_data=sqrt(N)*ifft(transmit_data_c2, N, 1);%ifft函数的操作按列进行，将时频域矩阵转换到时域.乘sqrt(N)保证功率归一化
        kron_c1 = kron(ones(1,Nsym),diag(Theta_c1'));
        transmit_data_noCP = kron_c1.*ifft_data;% 低复杂度计算
    end
    
    % 最后添加全CP
%     Ncp=floor(Tcp/T*N);%循环前缀的数据长度
    transmit_data_cp = zeros(Ncp,Nsym);
    for n=-Ncp:1:-1
        aa(Ncp+1+n) = exp(-1i*2*pi*c1*(N^2+2*N*n));
        transmit_data_cp(Ncp+1+n,:) = transmit_data_noCP(N+1+n,:).*exp(-1i*2*pi*c1*(N^2+2*N*n));
    end
    transmit_data_addCP=[transmit_data_cp;transmit_data_noCP];%，插入循环前缀，将每一帧数据的后端Ncp个数据补充到最前边
    % ifft_data1=[zeros(Ncp,Nsym);ifft_data];%，插入补零前缀，
    Tx_data=reshape(transmit_data_addCP,[],1);%并串转换，成为时域形式，发送   
end