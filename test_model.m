% 可选：初始化随机数生成器
% rng(0);

%% 生成测试数据
% rho_test = [100, 1000, 100, 10000, 2000]';
% z_test = [100, 2000, 4000, 8000, inf]';
% rho_test = [100, (10.^(2.5-1.5*square(2*pi.*linspace(0, 4, 48)))), 100]';
% z_test = [logspace(1, 4, 49), inf]';
% rho_test = (10.^(2.5-1.5*sin(linspace(0, 4*pi, 50))))';
% z_test = [logspace(1, 4, 49), inf]';
% rho_test = [10, 1000]';
% z_test = [1000, inf]';

% AQHK三层模型
rho_test = [10, 100, 1000]'; % A
% rho_test = [1000, 100, 10]'; % Q
% rho_test = [100, 10, 1000]'; % H
% rho_test = [100, 1000, 10]'; % K
z_test = [2000, 4000, inf]';

% 极端值测试模型
% rho_test = [250, 25, 100, 10, 25]';
% z_test = [600, 2000, 6000, 10000, inf]';

%% 正演
f_obs = logspace(-3, 3, 20)';
[rhoa_obs_log_clear, phs_obs_clear] = forward_func(log10(rho_test), log10(z_test), f_obs);

%% 生成误差
rhoa_obs_err_log = 0.1 * rhoa_obs_log_clear;
phs_obs_err = ones(length(f_obs), 1) * 5;
% rhoa_obs_err_log = rand(length(f_obs), 1) .* 0.3 .* rhoa_obs_log_clear;
% phs_obs_err = rand(length(f_obs), 1) .* 20;

%% 使相位与视电阻率数据量不对等
f_rhoa = f_obs;
f_phs = f_obs;
% f_phs = f_obs(1:2:end);
% phs_obs_clear = phs_obs_clear(1:2:end);
% phs_obs_err = phs_obs_err(1:2:end);

%% 生成噪声
% 均匀分布
% rhoa_obs_log = d_obs_log_clear + (2.*rand(length(f_obs), 1)-1).*rhoa_obs_err_log;
% phs_obs = phs_obs_clear + (2.*rand(length(f_obs), 1)-1).*phs_obs_err;
% 高斯误差
rhoa_obs_log = rhoa_obs_log_clear + randn(length(f_rhoa), 1).*(rhoa_obs_err_log);
phs_obs = phs_obs_clear + randn(length(f_phs), 1).*(phs_obs_err);
% 无
% rhoa_obs_log = rhoa_obs_log_clear;
% phs_obs = phs_obs_clear;