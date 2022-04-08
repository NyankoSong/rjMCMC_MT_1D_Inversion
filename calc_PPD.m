function model_grid = calc_PPD(model_cell, z_mesh, rho_mesh)
%CALC_PPD 计算后验概率密度
% 暂时未添加先验数据

model_len = length(model_cell); % 模型数
mesh_layers_n = length(z_mesh); % 网格层数
model_layers_vec_cell = cell(mesh_layers_n, 1);
lh_vec = [model_cell{:, 3}]'; % 似然度列向量
mesh_rho_n = length(rho_mesh); % 电阻网格数
model_grid = zeros(mesh_layers_n, mesh_rho_n); % 后验概率密度

% 重新整理成每层一个向量（TODO:低效，考虑删除）
for mesh_n_ind = 1:mesh_layers_n
    model_layers_vec_cell{mesh_n_ind} = zeros(mesh_layers_n, 1);
    for model_ind = 1:model_len
        zk = [model_cell{model_ind, 1}]-z_mesh(mesh_n_ind);
        z_ind = find(zk >= 0, 1);
        model_layers_vec_cell{mesh_n_ind}(model_ind) = model_cell{model_ind, 2}(z_ind);
    end
    disp(['整理中...', num2str(mesh_n_ind*100/mesh_layers_n), '%'])
end

% 按层计算后验概率密度（TODO:低效，考虑整合）
for mesh_n_ind = 1:mesh_layers_n
    model_layer_vec_k = model_layers_vec_cell{mesh_n_ind};
    model_layer_rho_k_unique = unique([model_layers_vec_cell{mesh_n_ind}])';
    for rho_k = model_layer_rho_k_unique
        model_grid(mesh_n_ind, find(rho_mesh == rho_k, 1)) = sum(lh_vec(model_layer_vec_k == rho_k));
    end
    disp(['计算中...', num2str(mesh_n_ind*100/mesh_layers_n), '%'])
end

% end

