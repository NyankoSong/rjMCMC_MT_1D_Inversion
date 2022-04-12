% 可选：初始化随机数生成器
% rng(0);

% 生成测试数据
m_test = [100, 1000, 100, 10000, 2000]';
z_test = [100, 2000, 4000, 8000, inf]';
% m_test = [100, (10.^(2.5-1.5*square(2*pi.*linspace(0, 4, 48)))), 100]';
% z_test = [logspace(1, 4, 49), inf]';
% m_test = (10.^(2.5-1.5*sin(linspace(0, 4*pi, 50))))';
% z_test = [logspace(1, 4, 49), inf]';
% m_test = [10, 1000]';
% z_test = [1000, inf]';

% 正演
f_obs = logspace(-3, 5, 20)';
[d_obs_log_clear, phs_obs_clear] = forward_func(log10(m_test), log10(z_test), f_obs);

% 生成误差
d_obs_err_log = 0.1 * d_obs_log_clear;
phs_obs_err = ones(length(f_obs), 1) * 5;

% 生成噪声
% 均匀分布
% d_obs_log = d_obs_log_clear + (2.*rand(length(f_obs), 1)-1).*d_obs_err_log;
% phs_obs = phs_obs_clear + (2.*rand(length(f_obs), 1)-1).*phs_obs_err;
% 正态分布
% d_obs_log = d_obs_log_clear + randn(length(f_obs), 1).*(d_obs_err_log/3);
% phs_obs = phs_obs_clear + randn(length(f_obs), 1).*(phs_obs_err/3);
% 无
d_obs_log = d_obs_log_clear;
phs_obs = phs_obs_clear;