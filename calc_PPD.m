function model_grid = calc_PPD(model_cell, z_mesh, rho_mesh)
%CALC_PPD �����������ܶ�
% ��ʱδ�����������

model_len = length(model_cell); % ģ����
mesh_layers_n = length(z_mesh); % �������
model_layers_vec_cell = cell(mesh_layers_n, 1);
lh_vec = [model_cell{:, 3}]'; % ��Ȼ��������
mesh_rho_n = length(rho_mesh); % ����������
model_grid = zeros(mesh_layers_n, mesh_rho_n); % ��������ܶ�

% ���������ÿ��һ��������TODO:��Ч������ɾ����
for mesh_n_ind = 1:mesh_layers_n
    model_layers_vec_cell{mesh_n_ind} = zeros(mesh_layers_n, 1);
    for model_ind = 1:model_len
        zk = [model_cell{model_ind, 1}]-z_mesh(mesh_n_ind);
        z_ind = find(zk >= 0, 1);
        model_layers_vec_cell{mesh_n_ind}(model_ind) = model_cell{model_ind, 2}(z_ind);
    end
    disp(['������...', num2str(mesh_n_ind*100/mesh_layers_n), '%'])
end

% ��������������ܶȣ�TODO:��Ч���������ϣ�
for mesh_n_ind = 1:mesh_layers_n
    model_layer_vec_k = model_layers_vec_cell{mesh_n_ind};
    model_layer_rho_k_unique = unique([model_layers_vec_cell{mesh_n_ind}])';
    for rho_k = model_layer_rho_k_unique
        model_grid(mesh_n_ind, find(rho_mesh == rho_k, 1)) = sum(lh_vec(model_layer_vec_k == rho_k));
    end
    disp(['������...', num2str(mesh_n_ind*100/mesh_layers_n), '%'])
end

% end

