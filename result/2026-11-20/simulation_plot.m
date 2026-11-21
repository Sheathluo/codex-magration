% close all
% clear

matFilePath = 'E:\Codex\codex-magration\result\2026-11-20\RMSE_result_20261120_181449.mat';
load(matFilePath);
% AFNE_delay_esti_RMSE(5) = 0.030434;
figure;
semilogy(para.SNR, AFNE_delay_esti_RMSE, 'b-o', 'LineWidth', 1.5);
hold on;
semilogy(para.SNR, RFT_delay_esti_RMSE, 'r+-', 'LineWidth', 1.5);
hold on;
semilogy(para.SNR, RFT_delay_esti_RMSE_old, 'g+-', 'LineWidth', 1.5);
grid on;
legend('AFNE','RFT','RFT-old');
xlabel('SNR (dB)');
ylabel('Normalized Delay RMSE');
title('Normalized Delay Estimation RMSE vs SNR');


% AFNE_doppler_esti_RMSE(5) = 0.0162;
figure;
semilogy(para.SNR, AFNE_doppler_esti_RMSE, 'b-o', 'LineWidth', 1.5);
hold on;
semilogy(para.SNR, RFT_doppler_esti_RMSE, 'r+-', 'LineWidth', 1.5);
hold on;
semilogy(para.SNR, RFT_doppler_esti_RMSE_old, 'g+-', 'LineWidth', 1.5);
grid on;
legend('AFNE','RFT','RFT-old');
xlabel('SNR (dB)');
ylabel('Normalized Doppler RMSE');
title('Normalized Doppler Estimation RMSE vs SNR');

