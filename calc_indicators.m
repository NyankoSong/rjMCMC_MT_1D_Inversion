model_grid_tmp = model_grid; % 选择要计算的model_grid

%% 计算峰值、期望值，统计层数
rho_average_log = sum(repmat(log10(rho_mesh)', length(z_mesh), 1) .* (model_grid_tmp ./ repmat(sum(model_grid_tmp, 2), 1, length(rho_mesh))), 2);
rho_average = 10.^rho_average_log;

[~, model_max_ind] = max(model_grid_tmp, [], 2);
rho_max = rho_mesh(model_max_ind);
rho_max_log = log10(rho_max);

model_n = [model_cell{:, 5}]';
model_n_hist = zeros(max(model_n) - min(model_n) + 1, 2);
for j = min(model_n):max(model_n)
    model_n_hist(j - min(model_n) + 1, :) = [j, sum(model_n(model_n == j)/j)];
end
model_n_hist(:, 2) = model_n_hist(:, 2) ./ sum(model_n_hist(:, 2));

%% 对峰值和期望值进行正演
[rhoa_average_log, ~] = forward_func(rho_average_log, z_mesh_log, f_rhoa);
[~, phs_average] = forward_func(rho_average_log, z_mesh_log, f_phs);
[rhoa_max_log, ~] = forward_func(rho_max_log, z_mesh_log, f_rhoa);
[~, phs_max] = forward_func(rho_max_log, z_mesh_log, f_phs);

%% 还原视电阻率
rhoa_average = 10.^rhoa_average_log;
rhoa_max = 10.^rhoa_max_log;

%% 计算模型中位数
model_grid_pdf = cumsum((model_grid_tmp ./ repmat(sum(model_grid_tmp, 2), 1, rho_n)), 2);
rho_median_log = zeros(z_n, 1);
for i = 1:z_n
    median_pos = find(model_grid_pdf(i, :) > 0.5, 1);
    rho_median_log(i) = interp1(model_grid_pdf(i, median_pos-1:median_pos), rho_mesh_log(median_pos-1:median_pos)', 0.5);
end
rho_median = 10.^rho_median_log;

%% 计算置信区间上界和下界
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

%% 计算模型标准差（用高斯分布描述置信区间时使用，然而实际上并不服从高斯分布......）
sigma = sum((model_grid_tmp ./ repmat(sum(model_grid_tmp, 2), 1, rho_n)) .* (repmat(rho_mesh_log', z_n, 1) - repmat(rho_average_log, 1, rho_n)).^2, 2) .^ (1/2);

%% 制图
plot_mesh