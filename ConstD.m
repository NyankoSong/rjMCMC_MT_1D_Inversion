function [model_cell, model_grid, end_flag, model_ind] = TransD(rho_mesh, z_mesh, f_obs, d_obs_log, d_obs_err_log, phase_obs, phase_obs_err, N, N_refresh, N_burn_in, rms_target, std_target, N_end, n)
% 单马尔科夫链程序——固定维

% N = 1E6; % 最大迭代次数
% N_refresh = 1E2; % 刷新间隔次数
% N_burn_in = 5E3; % 预热迭代次数
% rms_target = 0.2; % 目标RMS误差
% std_target = 0.1; % 期望标准差上限
% N_end = 5E3; % 判定终止范围

% 生成不确定性矩阵（TODO:整合进迭代，随迭代变化）
len = length(d_obs_log);
% Cd = (diag(linspace(1, 1, len), 0) + diag(linspace(0.5, 0.5, len-1), 1) + diag(linspace(0.5, 0.5, len-1), -1)) / 1E4;
Cd = diag(linspace(1, 1, len), 0) / 1E4;

% 处理网格
% z_mesh = logspace(0, 5, 200)';
% rho_mesh = logspace(0, 5, 200)';
z_mesh_log = log10(z_mesh);
rho_mesh_log = log10(rho_mesh);

rho_model_max = max(rho_mesh_log);

mesh_layers_n = length(z_mesh_log); % 网格层数
mesh_rho_n = length(rho_mesh_log); % 电阻网格数
model_grid = zeros(mesh_layers_n, mesh_rho_n); % 后验概率密度

% BOSITICK反演生成初始模型
[m_bostick, z_bostick] = bostick_func(d_obs_log, phase_obs, f_obs); % bostick反演
% n = 10; % 初始层数
n_range = [3, 30]; % 层数范围
% z = randi([1E0, 1E2], n, 1)*1E2;
z_log = [mesh_func(linspace(2, 4, n-1), z_mesh_log); inf];
% m = ceil(logspace(1, 4, n));
m_log = [interp1(z_bostick, m_bostick, z_log(1:end-1)); 0];
m_log(end) = m_log(end-1);
m_log = mesh_func(m_log, rho_mesh_log);
end_flag = 0;

cnt_rejections = 0;
model_ind = 1;
model_average_log_std = NaN;
% model_max_log = 0;
% model_average_log = 0;

lh = likelihood_func(Cd, d_obs_log, f_obs, m_log, z_log);
model_cell = cell(N, 4);
model_cell{1, 1} = z_log; % 分层向量
model_cell{1, 2} = m_log; % 电阻率向量
model_cell{1, 3} = lh; % 似然度
model_cell{1, 4} = cnt_rejections; % 拒绝次数
model_cell{1, 5} = n; % 层数
model_average_log_matrix = zeros(ceil(N/N_refresh), length(z_mesh));
refresh_ind = 0;

figure(1);
tic;
while model_ind < N
    flag = rand;
    m_log_new = m_log;
    z_log_new = z_log;
    n_new = n;
    if flag < 0.5 % 扰动电阻率 （TODO:考虑将操作概率转变为变量）
        n_rnd = randi(n_new);
        m_log_new(n_rnd) = proposal_func(m_log_new(n_rnd), 'rho', rho_mesh_log);
        if min(m_log_new) < rho_mesh_log(2)
            continue;
        end
    else % 扰动分层位置
        n_rnd = randi(n_new-1);
        z_log_new(n_rnd) = proposal_func(z_log_new(n_rnd), 'z', z_mesh_log);
        if z_log_new(1) < z_mesh_log(2) || min(z_log_new(2:end)-z_log_new(1:end-1)) < z_mesh_log(2)
            continue;
        end
%     else % 生灭
%         n_rnd = randi(n_new-2);
%         if (rand < 0.5 || n_new >= n_range(2)) && n_new > n_range(1) % 灭
%             z_log_new(n_rnd) = [];
%             m_log_new(n_rnd) = [];
%             n_new = n_new - 1;
%         else % 生
%             z_log_new = [z_log_new(1:n_rnd); mesh_func(mean(z_log_new(n_rnd:n_rnd+1)), z_mesh_log); z_log_new(n_rnd+1:end)];
%             m_log_new = [m_log_new(1:n_rnd); mesh_func(mean(m_log_new(n_rnd:n_rnd+1)), rho_mesh_log); m_log_new(n_rnd+1:end)];
%             n_new = n_new + 1;
%         end
    end
    lh_new = likelihood_func(Cd, d_obs_log, f_obs, m_log_new, z_log_new);
%     if rand < lh_new/lh && length(z_new) == length(unique(z_new)) && length(m_new) == length(unique(m_new)) && isequal(z_new, sort(z_new))
    if rand < lh_new/lh && length(z_log_new) == length(unique(z_log_new)) && isequal(z_log_new, sort(z_log_new)) && max(m_log_new) < rho_model_max
        model_ind = model_ind + 1;
        model_cell{model_ind, 1} = z_log_new;
        model_cell{model_ind, 2} = m_log_new;
        model_cell{model_ind, 3} = lh_new;
        model_cell{model_ind, 4} = cnt_rejections;
        model_cell{model_ind, 5} = n_new;
        m_log = m_log_new;
        z_log = z_log_new;
        lh = lh_new;
        n = n_new;
        cnt_rejections = 0;
        
        for mesh_n_ind = 1:mesh_layers_n % 计算后验概率
            zk = z_log - z_mesh_log(mesh_n_ind);
            z_ind = find(zk >= 0, 1);
            m_ind = find(rho_mesh_log == m_log(z_ind));
            model_grid(mesh_n_ind, m_ind) = model_grid(mesh_n_ind, m_ind) + lh;
        end
        
        if mod(model_ind, N_refresh) == 0
            refresh_ind = refresh_ind + 1;
            % 计算峰值和期望
            model_average_log = sum(repmat(rho_mesh_log', length(z_mesh_log), 1) .* (model_grid ./ repmat(sum(model_grid, 2), 1, length(rho_mesh_log))), 2);
            [~, model_max_log_ind] = max(model_grid, [], 2);
            model_max_log = rho_mesh_log(model_max_log_ind);
            
            % 判定期望的变化幅度
            model_average_log_matrix(refresh_ind, :) = model_average_log';
            if model_ind > N_end
                model_average_log_std = max(std(model_average_log_matrix(refresh_ind:-1:refresh_ind-round(N_end/N_refresh)+1, :)));
            end
                
            [rho_max_log, ~] = forward_func(model_max_log, z_mesh_log, f_obs);
            [rho_average_log, ~] = forward_func(model_average_log, z_mesh_log, f_obs);
            model_average = 10.^model_average_log;
            model_max = 10.^model_max_log;
            rms_max = rms((d_obs_log - rho_max_log)./d_obs_err_log);
            rms_average = rms((d_obs_log - rho_average_log)./d_obs_err_log);
            loglog(model_average, 10.^z_mesh_log)
            set(gca, 'YDir', 'reverse');
            axis([1, 1E5, 1, 1E5]);
            grid on
            hold on
            loglog(model_max, 10.^z_mesh_log, 'r')
            hold off
%             heatmap(model_grid)
%             grid off
            drawnow;
            disp(['迭代数：', num2str(model_ind), '，PPD峰值RMS误差：', num2str(rms_max), '，PPD期望RMS误差：', num2str(rms_average), '，std：', num2str(model_average_log_std)]);
            
            % 判断模型是否已稳定
            if model_average_log_std < std_target && model_ind > N_burn_in+N_end
                % 判断模型是否已达标
                if (rms_max < rms_target && rms_average < rms_target) || model_ind >= N
                    t = toc;
                    disp(['迭代结束，本次迭代', num2str(model_ind), '次，耗时', num2str(t), 's']);
                    break;
                end
                t = toc;
                disp(['模型已稳定，但不达标，迭代终止。本次迭代', num2str(model_ind), '次，耗时', num2str(t), 's']);
                end_flag = 1;
                break;
            end
        end
        
    else
        cnt_rejections = cnt_rejections + 1;
    end
end