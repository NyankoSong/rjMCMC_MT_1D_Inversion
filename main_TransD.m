%% 初始化参数
% 生成网格
z_n = 200;
rho_n = 200;

z_mesh = logspace(0, 5, z_n)';
rho_mesh = logspace(0, 6, rho_n)';
z_mesh_log = log10(z_mesh);
rho_mesh_log = log10(rho_mesh);

% 规定参数
N = 1E6; % 最大迭代次数
N_refresh = 1E3; % 刷新间隔次数
rms_target = 1; % 目标RMS误差
std_target = 0.05; % 期望标准差上限
N_end = 2E4; % 判定终止范围
k_punish = 1; % 罚参数（<1时为倾向更少层数）
f_Cd = 0.3; % 强相关频率差（暂定为1σ）
k_weight = 1; % 权重系数（<1时为视电阻率高权重）
k_err = 10; % 容差系数（越高则对误差容忍度越低）
z_smooth_log = 0.2; % 强相关层间距（暂定为1σ）

% 初始化参数
N_iter_sum = 0;
end_flag = 1;

n_test = 1;
t_mat = zeros(n_test, 1);
Niter_mat = zeros(n_test, 1);

%% 调用可变维函数
for i = 1:n_test
    
    t_main = tic;
    while end_flag == 1
        [model_cell, model_grid, end_flag, N_iter] = TransD(rho_mesh, z_mesh, f_obs, d_obs_log, d_obs_err_log, phs_obs, phs_obs_err, N, N_refresh, rms_target, std_target, N_end, k_punish, f_Cd, k_weight, k_err, z_smooth_log, m_test, z_test);
        N_iter_sum = N_iter_sum + N_iter;
        
        % 若只进行一次迭代则取消下行注释
        end_flag = 0;
    end
    t = toc(t_main);
    disp(['共迭代', num2str(N_iter_sum), '次，共耗时', num2str(t), 's'])
    
    end_flag = 1;
    t_mat(i) = t;
    Niter_mat(i) = N_iter_sum;
    N_iter_sum = 0;
end

%% 计算峰值、期望值，统计层数
model_average_log = sum(repmat(log10(rho_mesh)', length(z_mesh), 1) .* (model_grid ./ repmat(sum(model_grid, 2), 1, length(rho_mesh))), 2);
model_average = 10.^model_average_log;

[~, model_max_ind] = max(model_grid, [], 2);
model_max = rho_mesh(model_max_ind);
model_max_log = log10(model_max);

model_n = [model_cell{:, 5}]';
model_n_hist = zeros(max(model_n) - min(model_n) + 1, 2);
for j = min(model_n):max(model_n)
    model_n_hist(j - min(model_n) + 1, :) = [j, sum(model_n(model_n == j)/j)];
end
model_n_hist(:, 2) = model_n_hist(:, 2) ./ sum(model_n_hist(:, 2));

% 对峰值和期望值进行正演
[rhoa_average_log, phs_average] = forward_func(model_average_log, z_mesh_log, f_obs);
[rhoa_max_log, phs_max] = forward_func(model_max_log, z_mesh_log, f_obs);

% 还原视电阻率
rhoa_average = 10.^rhoa_average_log;
rhoa_max = 10.^rhoa_max_log;

% 计算模型中位数
model_grid_pdf = cumsum((model_grid ./ repmat(sum(model_grid, 2), 1, rho_n)), 2);
model_median_log = zeros(z_n, 1);
for i = 1:z_n
    median_pos = find(model_grid_pdf(i, :) > 0.5, 1);
    model_median_log(i) = interp1(model_grid_pdf(i, median_pos-1:median_pos), rho_mesh_log(median_pos-1:median_pos)', 0.5);
end
model_median = 10.^model_median_log;

% 计算置信区间上界和下界
model_conf = [0.99, 0.95, 0.68];
model_conf_edge_log = zeros(z_n, length(model_conf)*2);
for i = 1:z_n
    for j = 1:length(model_conf)
        conf_edge_pos_1 = find(model_grid_pdf(i, :) > (1-model_conf(j))/2, 1);
        conf_edge_pos_2 = find(model_grid_pdf(i, :) > (1-model_conf(j))/2 + model_conf(j), 1);
        if conf_edge_pos_1 > 1
            model_conf_edge_log(i, j*2-1) = interp1(model_grid_pdf(i, conf_edge_pos_1-1:conf_edge_pos_1), rho_mesh_log(conf_edge_pos_1-1:conf_edge_pos_1)', (1-model_conf(j))/2);
        else
            model_conf_edge_log(i, j*2-1) = rho_mesh_log(1);
        end
        model_conf_edge_log(i, j*2) = interp1(model_grid_pdf(i, conf_edge_pos_2-1:conf_edge_pos_2), rho_mesh_log(conf_edge_pos_2-1:conf_edge_pos_2)', (1-model_conf(j))/2 + model_conf(j));
    end
end
model_conf_edge = 10.^model_conf_edge_log;

% 计算模型标准差（用高斯分布描述置信区间时使用，然而实际上并不服从高斯分布......）
sigma = sum((model_grid ./ repmat(sum(model_grid, 2), 1, rho_n)) .* (repmat(rho_mesh_log', z_n, 1) - repmat(model_average_log, 1, rho_n)).^2, 2) .^ (1/2);

%% 制图
plot_mesh