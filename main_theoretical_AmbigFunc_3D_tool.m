%% 理论 pulse-shaped AFDM 模糊函数计算与三维显示工具
% 这个脚本独立于主蒙特卡洛链路，可用于：
%   1. 读取已有 MAT 缓存并画三维模糊函数；
%   2. 修改下面的参数后，重新计算理论模糊函数、保存缓存并画图。

clear;
close all;
% clc;
tic;

%% ===================== 0. 用户配置 =====================
cfg = struct();

% 计算/读取开关
cfg.recomputeAmbigFunc = true;     % true: 重新计算；false: 直接读取缓存
cfg.saveAmbigFunc = true;           % 重新计算后是否保存 MAT 缓存
cfg.ambigCachePath = './result/AmbiguFunc/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v6_20251022.mat';
cfg.ambigVarName = 'AmbigFunc_pulsedAFDM_noDelayDop';

% 图像显示控制
cfg.normalizeForPlot = true;        % 画图前是否按最大值归一化
cfg.dBFloor = -30;                  % dB 图动态范围下限
cfg.delayPlotStep = 1;              % delay 维抽点，1 表示不抽点
cfg.dopplerPlotStep = 1;            % Doppler 维抽点，1 表示不抽点

% 理论模糊函数参数。需要测试不同参数时，主要改这一组。
para = struct();
para.Nc = 128;
% para.Nc = 1024;
para.N = para.Nc;
para.L = 2;
para.beta = 0.35;
% para.beta = 0.25;
para.N_preSuffix = 20;
para.N_tapFiler = 20;
para.c1 = 8 / (2 * para.Nc);
% para.c1 = 11 / (2 * para.Nc);
para.mu_4 = 1.32;

cfg.L_delay = para.L;
para.L_Doppler = para.L;
% 理论模糊函数网格过采样倍率。
% delayAxis = 0 : 1/cfg.L_delay : para.N*para.L - 1/cfg.L_delay;
% dopplerAxis = -para.N*para.L/2 : 1/para.L_Doppler :para.N*para.L/2 - 1/para.L_Doppler;

% 理论模糊函数中心。当前主链路缓存使用零真实时延/多普勒。
cfg.tau0 = 0;
cfg.v0 = 0;

%% ===================== 1. 读取或重新计算模糊函数 =====================
if cfg.recomputeAmbigFunc
    fprintf('\nRecomputing theoretical ambiguity function...\n');
    fprintf('N = %d, L = %d, beta = %.3f, N_preSuffix = %d, N_tapFilter = %d\n' , ...
        para.N, para.L, para.beta, para.N_preSuffix, para.N_tapFiler);
    fprintf('rx_L_delay = %.6g, L_Doppler = %.6g\n', cfg.L_delay, para.L_Doppler);

    % 构造与主链路一致的等效 RRC 脉冲，并截断到一个 AFDM 符号周期内。
    % p_rrc_ambig = rcosdesign(para.beta, para.Nc, para.L, 'sqrt').';
    % p_rrc_ambig = p_rrc_ambig(1:end-1);

    % p_rrc_ambig_trunc = zeros(para.Nc * para.L, 1);
    % idx_start = para.Nc * para.L / 2 - para.N_preSuffix * para.L + 1;
    % idx_end = idx_start + (2 * para.N_preSuffix * para.L + 1) - 1;
    % p_rrc_ambig_trunc(idx_start:idx_end) = p_rrc_ambig(idx_start:idx_end);
    % p_rrc_ambig_trunc = p_rrc_ambig
    % p_rrc_ambig = circshift(p_rrc_ambig_trunc, -para.Nc * para.L / 2);

    % 为了与使用中的脉冲成形滤波器保持一致
    p_rrc = rcosdesign(para.beta, para.N_tapFiler, para.L, 'sqrt').';
    p_rrc_trunc = zeros(para.Nc * para.L, 1);
    start_idx_rrc = para.Nc * para.L / 2 - para.N_preSuffix * para.L / 2 + 1;
    end_idx_rrc = start_idx_rrc + length(p_rrc) - 1;

    p_rrc_trunc(start_idx_rrc:end_idx_rrc) = p_rrc;
    p_rrc_ambig = circshift(p_rrc_trunc, -para.Nc * para.L / 2);

    delayAxis = 0:1/cfg.L_delay:para.N * para.L - 1/cfg.L_delay;
    dopplerAxis = -para.N * para.L / 2:1/para.L_Doppler:para.N * para.L / 2 - 1/para.L_Doppler;

    % 输出变量名与主链路保持一致，方便后续直接复用缓存。
    AmbigFunc_pulsedAFDM_noDelayDop = f_theoretical_AmbigFunc_pulseShaped_AFDM_v7( ...
        p_rrc_ambig, para.N, para.L, para.c1, para.mu_4, delayAxis, dopplerAxis, cfg.tau0, cfg.v0);

    if cfg.saveAmbigFunc
        cacheDir = fileparts(cfg.ambigCachePath);
        if ~isempty(cacheDir) && ~exist(cacheDir, 'dir')
            mkdir(cacheDir);
        end

        % 重新计算得到的新模糊函数不覆盖原缓存文件，而是按当前日期时间另存一份
        [cacheDir, cacheBaseName, cacheExt] = fileparts(cfg.ambigCachePath);
        timeStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        ambigSavePath = fullfile(cacheDir, [cacheBaseName, '_', timeStamp, cacheExt]);

        save(ambigSavePath, 'AmbigFunc_pulsedAFDM_noDelayDop', 'para', 'cfg', ...
            'delayAxis', 'dopplerAxis');
        fprintf('Saved ambiguity function to: %s\n', ambigSavePath);
    end
else
    if ~exist(cfg.ambigCachePath, 'file')
        error('未找到模糊函数缓存文件: %s', cfg.ambigCachePath);
    end

    cacheInfo = whos('-file', cfg.ambigCachePath);
    cacheVarNames = string({cacheInfo.name});
    if ~any(cacheVarNames == cfg.ambigVarName)
        fprintf('缓存文件中包含的变量如下：\n');
        disp(cacheVarNames.');
        error('缓存文件中没有变量 %s，请检查 cfg.ambigVarName。', cfg.ambigVarName);
    end

    ambigData = load(cfg.ambigCachePath, cfg.ambigVarName);
    AmbigFunc_pulsedAFDM_noDelayDop = ambigData.(cfg.ambigVarName);

    % 优先使用缓存里的坐标轴；没有坐标轴时再根据矩阵尺寸和当前 N/L 反推。
    cacheAxisInfo = whos('-file', cfg.ambigCachePath, 'delayAxis', 'dopplerAxis');
    if numel(cacheAxisInfo) == 2
        axisData = load(cfg.ambigCachePath, 'delayAxis', 'dopplerAxis');
        delayAxis = axisData.delayAxis;
        dopplerAxis = axisData.dopplerAxis;
    else
        delayAxis = [];
        dopplerAxis = [];
    end
end

% ambigCachePath_real = './result/AmbiguFunc/AmbigFunc_pulsedAFDM_th_fracDelayFracDop_N128_2Nc1_8_L2_Ldelay2_Ldoppler2_v6_20251022.mat';
% ambigVarName_real = 'AmbigFunc_pulsedAFDM_noDelayDop';
% ambigData_real = load(ambigCachePath_real, ambigVarName_real);
% AmbigFunc_pulsedAFDM_noDelayDop_real = ambigData_real.(ambigVarName_real);
% err = sum(abs(AmbigFunc_pulsedAFDM_noDelayDop_real-AmbigFunc_pulsedAFDM_noDelayDop),"all")

%% ===================== 2. 坐标轴与幅度处理 =====================
AmbigAbs = abs(AmbigFunc_pulsedAFDM_noDelayDop);
if cfg.normalizeForPlot
    AmbigPlotBase = AmbigAbs ./ (max(AmbigAbs(:)) + eps);
    zLabelLinear = 'Normalized power';
else
    AmbigPlotBase = AmbigAbs;
    zLabelLinear = 'Power';
end
% 理论函数输出 E{|A(tau,nu)|^2}，属于功率量，因此使用 10*log10。
AmbigPlotBase_dB = 10 * log10(AmbigPlotBase ./ (max(AmbigPlotBase(:)) + eps) + eps);

[numDelay, numDoppler] = size(AmbigPlotBase);

if isempty(delayAxis) || isempty(dopplerAxis)
    delayOversample = numDelay / (para.N * para.L);
    dopplerOversample = numDoppler / (para.N * para.L);
    if abs(delayOversample - round(delayOversample)) > 1e-12 || ...
            abs(dopplerOversample - round(dopplerOversample)) > 1e-12
        warning('矩阵尺寸无法由 N*L 整除，将使用索引坐标画图。');
        delayAxis = 0:numDelay-1;
        dopplerAxis = 0:numDoppler-1;
        delayLabel = 'Delay index';
        dopplerLabel = 'Doppler index';
    else
        delayAxis = 0:1/delayOversample:para.N * para.L - 1/delayOversample;
        dopplerAxis = -para.N * para.L / 2:1/dopplerOversample:para.N * para.L / 2 - 1/dopplerOversample;
        delayLabel = 'Normalized delay';
        dopplerLabel = 'Normalized Doppler';
    end
else
    delayOversample = numDelay / (para.N * para.L);
    dopplerOversample = numDoppler / (para.N * para.L);
    delayLabel = 'Normalized delay';
    dopplerLabel = 'Normalized Doppler';
end

if numel(delayAxis) ~= numDelay || numel(dopplerAxis) ~= numDoppler
    error('坐标轴长度与模糊函数矩阵尺寸不一致。delayAxis=%d, dopplerAxis=%d, matrix=%dx%d', ...
        numel(delayAxis), numel(dopplerAxis), numDelay, numDoppler);
end

delayIdx = 1:cfg.delayPlotStep:numDelay;
dopplerIdx = 1:cfg.dopplerPlotStep:numDoppler;
delayAxisPlot = delayAxis(delayIdx);
dopplerAxisPlot = dopplerAxis(dopplerIdx);
AmbigPlot = AmbigPlotBase(delayIdx, dopplerIdx);
AmbigPlot_dB = AmbigPlotBase_dB(delayIdx, dopplerIdx);

% 时延维是长度为 N*L 的循环坐标。fftshift 后用有符号时延轴显示，
% 将原来分布在首尾两端的零时延主瓣重新拼接到图中央。
if numDelay > 1
    delayGridSpacing = median(diff(delayAxis));
else
    delayGridSpacing = 1;
end
delayAxisShift = (-floor(numDelay / 2):ceil(numDelay / 2)-1) * delayGridSpacing;
AmbigPlotBaseShift = fftshift(AmbigPlotBase, 1);
AmbigPlotBaseShift_dB = fftshift(AmbigPlotBase_dB, 1);

delayShiftIdx = 1:cfg.delayPlotStep:numDelay;
delayAxisShiftPlot = delayAxisShift(delayShiftIdx);
AmbigPlotShift = AmbigPlotBaseShift(delayShiftIdx, dopplerIdx);
AmbigPlotShift_dB = AmbigPlotBaseShift_dB(delayShiftIdx, dopplerIdx);

fprintf('\n=========== Ambiguity Function Info ===========\n');
fprintf('Mode       : %s\n', ternaryString(cfg.recomputeAmbigFunc, 'recompute', 'load cache'));
fprintf('File       : %s\n', cfg.ambigCachePath);
fprintf('Matrix size: %d delay bins x %d Doppler bins\n', numDelay, numDoppler);
fprintf('Peak value : %.6g before plot normalization\n', max(AmbigAbs(:)));
fprintf('Delay oversampling   : %.6g\n', delayOversample);
fprintf('Doppler oversampling : %.6g\n', dopplerOversample);
fprintf('Elapsed time: %.2f s\n', toc);
fprintf('================================================\n\n');

%% ===================== 3. 未移位的三维图与俯视图 =====================
[DopplerGrid, DelayGrid] = meshgrid(dopplerAxisPlot, delayAxisPlot);

figure('Name', 'Ambiguity Function - Linear Scale - Unshifted Delay', 'Color', 'w');
surf(DopplerGrid, DelayGrid, AmbigPlot, 'EdgeColor', 'none');
shading interp;
colormap turbo;
colorbar;
xlabel(dopplerLabel);
ylabel(delayLabel);
zlabel(zLabelLinear);
title('Pulse-shaped AFDM Ambiguity Function Power |A(\tau,\nu)|^2 (Unshifted Delay)');
view(45, 35);
grid on;

figure('Name', 'Ambiguity Function - dB Scale - Unshifted Delay', 'Color', 'w');
surf(DopplerGrid, DelayGrid, AmbigPlot_dB, 'EdgeColor', 'none');
shading interp;
colormap turbo;
colorbar;
xlabel(dopplerLabel);
ylabel(delayLabel);
zlabel('Power (dB)');
title('Pulse-shaped AFDM Ambiguity Function Power in dB (Unshifted Delay)');
view(45, 35);
grid on;
zlim([cfg.dBFloor, 0]);

figure('Name', 'Ambiguity Function - Top View', 'Color', 'w');
imagesc(dopplerAxisPlot, delayAxisPlot, AmbigPlot_dB);
axis xy;
colormap turbo;
colorbar;
clim([cfg.dBFloor, 0]);
xlabel(dopplerLabel);
ylabel(delayLabel);
title('Pulse-shaped AFDM Ambiguity Function Top View (dB, Unshifted Delay)');

%% ===================== 4. 时延维 fftshift 后的三维图与俯视图 =====================
[DopplerGridShift, DelayGridShift] = meshgrid(dopplerAxisPlot, delayAxisShiftPlot);

figure('Name', 'Ambiguity Function - Linear Scale - Shifted Delay', 'Color', 'w');
surf(DopplerGridShift, DelayGridShift, AmbigPlotShift, 'EdgeColor', 'none');
shading interp;
colormap turbo;
colorbar;
xlabel(dopplerLabel);
ylabel([delayLabel, ' (centered)']);
zlabel(zLabelLinear);
title('Pulse-shaped AFDM Ambiguity Function Power |A(\tau,\nu)|^2 (fftshift on Delay)');
view(45, 35);
grid on;

figure('Name', 'Ambiguity Function - dB Scale - Shifted Delay', 'Color', 'w');
surf(DopplerGridShift, DelayGridShift, AmbigPlotShift_dB, 'EdgeColor', 'none');
shading interp;
colormap turbo;
colorbar;
xlabel(dopplerLabel);
ylabel([delayLabel, ' (centered)']);
zlabel('Power (dB)');
title('Pulse-shaped AFDM Ambiguity Function in dB (fftshift on Delay)');
view(45, 35);
grid on;
zlim([cfg.dBFloor, 0]);

figure('Name', 'Ambiguity Function - Top View - Shifted Delay', 'Color', 'w');
imagesc(dopplerAxisPlot, delayAxisShiftPlot, AmbigPlotShift_dB);
axis xy;
colormap turbo;
colorbar;
clim([cfg.dBFloor, 0]);
xlabel(dopplerLabel);
ylabel([delayLabel, ' (centered)']);
title('Pulse-shaped AFDM Ambiguity Function Top View (dB, fftshift on Delay)');

%% ===================== 5. 未移位与移位后的切片 =====================
% 分别在理论模糊函数中心多普勒和中心时延处取最近网格点，
% 绘制穿过主峰的时延切片与多普勒切片。
[~, dopplerSliceIdx] = min(abs(dopplerAxis - cfg.v0));
[~, delaySliceIdx] = min(abs(delayAxis - cfg.tau0));

delaySlice_dB = AmbigPlotBase_dB(delayIdx, dopplerSliceIdx);
delaySliceShift_dB = AmbigPlotBaseShift_dB(delayShiftIdx, dopplerSliceIdx);
dopplerSlice_dB = AmbigPlotBase_dB(delaySliceIdx, dopplerIdx);

figure('Name', 'Ambiguity Function - Delay Slice', 'Color', 'w');
plot(delayAxisPlot, delaySlice_dB, 'LineWidth', 1.5);
grid on;
xlabel(delayLabel);
ylabel('Power (dB)');
title(sprintf('Pulse-shaped AFDM Delay Slice (Doppler = %.6g, Unshifted)', ...
    dopplerAxis(dopplerSliceIdx)));
ylim([cfg.dBFloor, 0]);

figure('Name', 'Ambiguity Function - Delay Slice - Shifted', 'Color', 'w');
plot(delayAxisShiftPlot, delaySliceShift_dB, 'LineWidth', 1.5);
grid on;
xlabel([delayLabel, ' (centered)']);
ylabel('Power (dB)');
title(sprintf('Pulse-shaped AFDM Delay Slice (Doppler = %.6g, fftshift)', ...
    dopplerAxis(dopplerSliceIdx)));
ylim([cfg.dBFloor, 0]);

figure('Name', 'Ambiguity Function - Doppler Slice', 'Color', 'w');
plot(dopplerAxisPlot, dopplerSlice_dB, 'LineWidth', 1.5);
grid on;
xlabel(dopplerLabel);
ylabel('Power (dB)');
title(sprintf('Pulse-shaped AFDM Doppler Slice (Delay = %.6g)', ...
    delayAxis(delaySliceIdx)));
ylim([cfg.dBFloor, 0]);

%% ===================== 本脚本局部辅助函数 =====================
function out = ternaryString(cond, trueStr, falseStr)
    if cond
        out = trueStr;
    else
        out = falseStr;
    end
end
