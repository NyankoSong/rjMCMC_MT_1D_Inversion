%% 初始化参数
% 生成网格
z_n = 100;
rho_n = 50;

z_mesh = logspace(0, 4, z_n)';
rho_mesh = logspace(0, 4, rho_n)';
z_mesh_log = log10(z_mesh);
rho_mesh_log = log10(rho_mesh);

% 规定参数
N = 1E6; % 最大迭代次数
N_refresh = 1E3; % 刷新间隔次数（同时也是误差标准差的采样长度）
rms_target = 1; % 目标RMS误差
burn_in_std_target = 1; % 预热阶段误差标准差的标准差上限
std_target = 0.01; % 误差标准差的标准差上限
N_end = 1E4; % 判定终止范围（必须是N_refresh的整数倍）
k_punish = 1; % 罚参数（<1时为倾向更少层数）
k_weight = 1; % 权重系数（<1时为视电阻率高权重）
burn_in_data_scale_factor = 2; % 预热阶段数据协方差矩阵尺度因子
main_data_scale_factor = 1; % 数据协方差矩阵尺度因子
model_scale_factor = 1; % 模型协方差矩阵尺度因子
lambda = 1; % 先验权系数
z_smooth_log = 0; % 强相关层间距（暂定为1σ）

%% 调用可变维函数
t_main = tic;
[model_cell, burn_in_model_grid, model_grid, end_flag, N_iter, N_iter_burn_in] = TransD(rho_mesh, z_mesh, f_rhoa, f_phs, rhoa_obs_log, rhoa_obs_err_log, phs_obs, phs_obs_err, N, N_refresh, rms_target, burn_in_std_target, std_target, N_end, k_punish, k_weight, burn_in_data_scale_factor, main_data_scale_factor, model_scale_factor, lambda, z_smooth_log, rho_test, z_test);
t = toc(t_main);

if end_flag == 0
%     disp(['模型达标，共迭代', num2str(N_iter), '次，共耗时', num2str(t), 's'])
elseif end_flag == 1
%     disp(['模型不达标，共迭代', num2str(N_iter), '次，共耗时', num2str(t), 's'])
elseif end_flag == 2
    disp('已拒绝模型数达到上限，终止运行')
end

%% 计算各项指标
calc_indicators