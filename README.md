# AFDM 距离徙动补偿与目标参数估计

本项目使用 MATLAB 仿真脉冲成形的仿射频分复用（AFDM）通感一体化波形，研究目标在相干积累期间发生距离徙动时的回波补偿，以及时延、多普勒和回波增益估计。仿真流程包括 AFDM 发射信号生成与根升余弦（RRC）成形、时变信道回波生成、接收预处理、匹配滤波粗搜索和牛顿迭代精估计。仓库也提供 Radon–Fourier 变换（RFT）实现和已保存的仿真结果。

## 主要入口

| 文件 | 用途 |
| --- | --- |
| [`main_20261104_v11_MagriCom_AFDM_N128.m`](main_20261104_v11_MagriCom_AFDM_N128.m) | 单次回波仿真，运行多目标 AFNE 粗搜索与牛顿正交匹配追踪（OMP）精估计。 |
| [`main_20261114_v12_single_target_Monte_MagriCom_N128.m`](main_20261114_v12_single_target_Monte_MagriCom_N128.m) | 单目标 AFNE 蒙特卡洛仿真，按信噪比统计归一化时延、多普勒和回波增益误差。 |
| [`main_theoretical_AmbigFunc_3D_tool.m`](main_theoretical_AmbigFunc_3D_tool.m) | 独立计算或读取理论模糊函数，并绘制三维曲面、俯视图和切片。 |
| [`RFT.m`](RFT.m)、[`RFT_v2.m`](RFT_v2.m)、[`f_AFDM_RFT_v5.m`](f_AFDM_RFT_v5.m) | RFT 相关函数。`RFT_v2.m` 在脉压数据的距离维进行分数格取样，`f_AFDM_RFT_v5.m` 对 AFDM 接收符号进行徙动补偿和相干积累。 |
| [`result/`](result/) | 已保存的 MAT 仿真结果、绘图脚本、参数记录及部分接收采样数据。 |

## 运行仿真

1. 在 MATLAB 中打开项目根目录。主链路需要 `qammod`、`bi2de` 和 `rcosdesign` 等函数，可用 `which` 检查本机的安装情况。
2. 准备与当前波形和搜索网格匹配的理论模糊函数 MAT 文件。主脚本通过 `cachePath` 指定它，并读取变量 `AmbigFunc_pulsedAFDM_noDelayDop`。`AmbiguFunc/` 是本地缓存目录，克隆后需要自行准备。
3. 如需生成缓存，先在 `main_theoretical_AmbigFunc_3D_tool.m` 中启用 `cfg.recomputeAmbigFunc` 和 `cfg.saveAmbigFunc`，核对 AFDM 符号长度、啁啾参数、RRC 脉冲成形与时延和多普勒采样网格，并将 `cfg.ambigCachePath` 指向本地 `AmbiguFunc/`。工具会将缓存另存为带时间戳的 MAT 文件。将实际生成的路径填入主脚本的 `cachePath`。
4. 运行单次仿真，确认参数和搜索范围后，再运行蒙特卡洛脚本。可在脚本开头设置 `para.SNR`、`para.Ntrial` 和 `para.Nsym`，在感知配置中设置目标参数范围和 `cfg.flag_delayMigration`（`1` 为符号级徙动，`2` 为采样点级徙动，`3` 为关闭徙动）。

```matlab
% 在项目根目录运行。先在脚本中设置正确的 cachePath。
main_20261104_v11_MagriCom_AFDM_N128
```

部分历史脚本依赖旧版模糊函数读取函数或本地硬件数据。新环境下可从上表的主脚本和已保存结果开始。

## 查看已保存的结果

蒙特卡洛脚本将配置、真实参数、估计值和均方根误差（RMSE）保存到 `result/年-月-日/RMSE_result_时间戳.mat`。在项目根目录可直接读取已有结果并绘图，无需生成模糊函数缓存。

```matlab
data = load(fullfile('result', '2026-11-16', 'RMSE_result_20261116_212424.mat'));
semilogy(data.para.SNR, data.AFNE_delay_esti_RMSE, '-o');
grid on;
xlabel('信噪比 (dB)');
ylabel('归一化时延 RMSE');
```

当前单目标蒙特卡洛脚本的 RFT 调用处于注释状态，相关输出使用占位值。做 AFNE 与 RFT 的性能对比时，需要启用实际 RFT 估计，并统一波形、信噪比、目标参数范围和误差单位。
