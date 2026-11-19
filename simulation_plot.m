% close all
% clear

sample_data_AFNE = load("result/2026-09-24/AFDM_MC_SNR_RMSE_N128_20260924_144322.mat");
symbol_data_RFT = load("result/2026-09-24/AFDM_MC_SNR_RMSE_N128_20260924_144322.mat");
symbol_sim_result_AFNE = sample_data_AFNE.mc_result;
RFT_symbol_sim_result = symbol_data_RFT.mc_result;

figure;
semilogy(RFT_symbol_sim_result.SNR_dB_set, RFT_symbol_sim_result.RFT_delay_rmse, 'b-o', 'LineWidth', 1.5);
hold on;
semilogy(symbol_sim_result_AFNE.SNR_dB_set, symbol_sim_result_AFNE.AFNE_delay_rmse, 'r+-', 'LineWidth', 1.5);
grid on;
legend('RFT','AFNE');
legend('AFNE');
xlabel('SNR (dB)');
ylabel('Normalized Delay RMSE');
title('Normalized Delay Estimation RMSE vs SNR');

figure;
semilogy(RFT_symbol_sim_result.SNR_dB_set, RFT_symbol_sim_result.RFT_doppler_rmse, 'b-o', 'LineWidth', 1.5);
hold on;
semilogy(symbol_sim_result_AFNE.SNR_dB_set, symbol_sim_result_AFNE.AFNE_doppler_rmse, 'r+-', 'LineWidth', 1.5);
grid on;
% legend('RFT','AFNE');
legend('AFNE');
xlabel('SNR (dB)');
ylabel('Normalized Doppler RMSE');
title('Normalized Doppler Estimation RMSE vs SNR');

% figure;
% semilogy(RFT_symbol_sim_result.SNR_dB_set, RFT_symbol_sim_result.RFT_range_rmse_m, 'b-o', 'LineWidth', 1.5);
% hold on;
% semilogy(symbol_sim_result_AFNE.SNR_dB_set, symbol_sim_result_AFNE.AFNE_range_rmse_m, 'r+-', 'LineWidth', 1.5);
% grid on;
% legend('RFT','AFNE');
% xlabel('SNR (dB)');
% ylabel('Range RMSE (m)');
% title('Range Estimation RMSE vs SNR');
%
% figure;
% semilogy(RFT_symbol_sim_result.SNR_dB_set, RFT_symbol_sim_result.RFT_velocity_rmse_mps, 'b-o', 'LineWidth', 1.5);
% hold on;
% semilogy(symbol_sim_result_AFNE.SNR_dB_set, symbol_sim_result_AFNE.AFNE_velocity_rmse_mps, 'r+-', 'LineWidth', 1.5);
% grid on;
% legend('RFT','AFNE');
% xlabel('SNR (dB)');
% ylabel('Velocity RMSE (m/s)');
% title('Velocity Estimation RMSE vs SNR');
