%% Starlink-like pass over Beijing: relative velocity with ground station
% 轨道: 70 deg / 570 km (Starlink 高纬壳, 可天顶过北京)
% 物理: 圆轨道 + 地球自转(WGS84地面站) + ECI->ECEF 变换
%       相对速度 = v_sat(ECEF) - v_station;  多普勒分量 = 视线径向
% 用法: 直接运行, 生成 3 张独立 figure

clear; clc; close all;

%% ===== 1. 物理常数 =====
mu      = 398600.4418;     % 地球引力常数 km^3/s^2
Rearth  = 6378.137;        % WGS84 赤道半径 km
f_flat  = 1/298.257223563; % WGS84 扁率
e2      = f_flat*(2-f_flat);
omega_e = 7.292115e-5;     % 地球自转角速度 rad/s
c_light = 299792.458;      % 光速 km/s

%% ===== 2. 轨道参数 (圆轨道) =====
h      = 570;              % 轨道高度 km
a      = Rearth + h;       % 半长轴 km
i_deg  = 70;               % 轨道倾角 deg
i      = deg2rad(i_deg);
w      = 0;                % 近地点幅角 (圆轨道可设 0)

%% ===== 3. 北京地面站 (WGS84) =====
latB = deg2rad(39.9042);
lonB = deg2rad(116.4074);
altB = 0.05;               % 海拔 km
NB   = Rearth/sqrt(1-e2*sin(latB)^2);
rB_ecef = [ (NB+altB)*cos(latB)*cos(lonB);
            (NB+altB)*cos(latB)*sin(lonB);
            (NB*(1-e2)+altB)*sin(latB) ];      % 站址 ECEF (km)
omega_vec = [0;0;omega_e];
vB_ecef   = cross(omega_vec, rB_ecef);         % 地面站速度 (km/s)

%% ===== 4. 求解"天顶过北京"的初始条件 (RAAN, nu) =====
% 通过数值最小化, 使 t=0 时卫星星下点 = 北京经纬度
Ro = @(Om) [ cos(Om) -sin(Om) 0; sin(Om) cos(Om) 0; 0 0 1];
Ri = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
Rw = [cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1];

obj = @(x) latlon_residual(x, a, Ri, Rw, Ro, latB, lonB);
x0  = [deg2rad(90); deg2rad(45)];
opts = optimset('Display','off','TolX',1e-9,'TolFun',1e-15);
xopt  = fminsearch(obj, x0, opts);
Om_val = xopt(1);  nu_center = xopt(2);
R = Ro(Om_val)*Ri*Rw;   % 总旋转矩阵 (轨道面 -> ECI)

%% ===== 5. 时间轴 (以天顶过境为中心 ±300 s) =====
n_orbit = sqrt(mu/a^3);          % 平运动角速度 rad/s
T       = 2*pi/n_orbit;          % 周期 s
tspan   = -300:1:300;
nu_ref  = nu_center + tspan*n_orbit;

v0 = sqrt(mu/a);                 % 圆轨道速度大小 km/s

%% ===== 6. 主循环: 逐时刻 ECI->ECEF, 算相对速度 =====
Nt = numel(tspan);
vrel_mag = zeros(1,Nt);
vrad     = zeros(1,Nt);
elv      = zeros(1,Nt);
rng_km   = zeros(1,Nt);

for k = 1:Nt
    nu = nu_ref(k);
    r_orb = a*[cos(nu); sin(nu); 0];
    v_orb = v0*[-sin(nu); cos(nu); 0];

    r_eci = R*r_orb;
    v_eci = R*v_orb;

    theta = omega_e*tspan(k);                 % 简化恒星时 (相对演示)
    Rz = [cos(theta) sin(theta) 0;            % ECI -> ECEF  (IERS 惯例)
         -sin(theta) cos(theta) 0;
          0           0         1];

    r_ecef = Rz*r_eci;
    v_ecef = Rz*v_eci - cross(omega_vec, r_ecef);

    rho   = r_ecef - rB_ecef;
    rng_km(k) = norm(rho);

    vrel      = v_ecef - vB_ecef;
    vrel_mag(k) = norm(vrel);
    vrad(k)     = dot(vrel, rho)/norm(rho);   % 视线径向分量

    [S,E,Z] = ecef_to_sez(rho, latB, lonB);
    elv(k)  = atan2(Z, sqrt(S^2+E^2))*180/pi;
end

%% ===== 7. 输出关键结果 =====
[~,idx] = max(elv);
fprintf('===== 轨道: i=%g deg, h=%g km, a=%g km =====\n', i_deg, h, a);
fprintf('周期 T = %.1f s = %.2f min\n', T, T/60);
fprintf('轨道速度 v0 = %.0f m/s,  地面站自转速度 = %.0f m/s\n', ...
        v0*1000, norm(vB_ecef)*1000);
fprintf('天顶时刻 t = %g s,  仰角 = %.2f deg (目标 90)\n', tspan(idx), elv(idx));
fprintf('  相对速率 |v_rel| = %.0f m/s\n', vrel_mag(idx)*1000);
fprintf('  径向速度(多普勒) = %.0f m/s  (天顶应接近0)\n', vrad(idx)*1000);
fprintf('  斜距 = %.0f km (正天顶应≈h=%g)\n', rng_km(idx), h);

fc = 12e9;   % 示例载波 Ku/Ka 段
fprintf('示例 fc=%.0f GHz -> 多普勒 = %.1f kHz\n', fc/1e9, vrad(idx)/c_light*fc/1e3);

%% ===== 8. 绘图: 三张独立 figure =====
% --- 图1: 仰角 ---
figure('Name','北京站仰角','NumberTitle','off','Position',[100 100 640 400]);
plot(tspan, elv, 'b-', 'LineWidth', 1.3); hold on;
plot(tspan(idx), elv(idx), 'ro', 'MarkerSize', 9, 'LineWidth', 1.5);
grid on; xlabel('t 相对天顶过境 (s)'); ylabel('仰角 (deg)');
title(sprintf('北京站仰角   (i=%g°, h=%gkm)', i_deg, h));
legend('仰角曲线','天顶时刻','Location','North');

% --- 图2: 相对速率 ---
figure('Name','相对速率','NumberTitle','off','Position',[780 100 640 400]);
plot(tspan, vrel_mag*1000, 'k-', 'LineWidth', 1.3); hold on;
plot(tspan(idx), vrel_mag(idx)*1000, 'ro', 'MarkerSize', 9, 'LineWidth', 1.5);
grid on; xlabel('t (s)'); ylabel('相对速率 (m/s)');
title('卫星—接收端相对速率 |v_{rel}|');
legend('相对速率','天顶时刻','Location','North');

% --- 图3: 视线径向速度 / 多普勒 ---
figure('Name','多普勒径向速度','NumberTitle','off','Position',[100 560 640 400]);
plot(tspan, vrad*1000, 'm-', 'LineWidth', 1.3); hold on;
plot(tspan(idx), vrad(idx)*1000, 'ro', 'MarkerSize', 9, 'LineWidth', 1.5);
grid on; xlabel('t (s)'); ylabel('径向速度 (m/s)');
title('视线径向速度 v_{rad} (多普勒分量, 穿越零点)');
legend('径向速度','天顶时刻','Location','South');

%% ===== 子函数 =====
function f = latlon_residual(x, a, Ri, Rw, Ro, latB, lonB)
    % 残差: 星下点经纬度 与 北京 的差异
    Om = x(1); nu = x(2);
    R  = Ro(Om)*Ri*Rw;
    r  = a*[cos(nu); sin(nu); 0];
    r_eci = R*r;
    lat = asin(r_eci(3)/a);
    lon = atan2(r_eci(2), r_eci(1));
    dlon = atan2(sin(lon-lonB), cos(lon-lonB));   % 最短经度差
    f = (lat-latB)^2 + dlon^2;
end

function [S,E,Z] = ecef_to_sez(rho, lat, lon)
    % ECEF -> 站心 SEZ (South-East-Zenith)
    Rse = [ sin(lat)*cos(lon)  sin(lat)*sin(lon) -cos(lat);
           -sin(lon)           cos(lon)           0;
            cos(lat)*cos(lon)  cos(lat)*sin(lon)  sin(lat)];
    sez = Rse*rho;
    S = sez(1); E = sez(2); Z = sez(3);
end
