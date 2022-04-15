function [model_cell, model_grid, end_flag, model_ind] = TransD(rho_mesh, z_mesh, f_rhoa, f_phs, rhoa_obs_log, rhoa_obs_err_log, phs_obs, phs_obs_err, N, N_refresh, rms_target, std_target, N_end, k_punish, k_weight, scale_factor, lambda, z_smooth_log, m_test, z_test)
% 单马尔科夫链程序——可变维
% end_flag % 迭代结束时是否达标

% N = 1E6; % 最大迭代次数
% N_refresh = 1E2; % 刷新间隔次数
% N_burn_in = 5E3; % 预热迭代次数
% rms_target = 0.2; % 目标RMS误差
% std_target = 0.1; % 期望标准差上限
% N_end = 5E3; % 判定终止范围
% k_punish = 0.9; % 罚参数
% k_weight = 1; % 权重系数（<1时为视电阻率高权重）
% k_err = 1E4; % 容差系数
% z_smooth_log = 0.2; % 强相关层间距（暂定为1σ）

print_flag = 0; % 是否输出帧到./Frames/tmp
plot_flag = 1; % 是否实时进行制图

if length(m_test) > 1
    test_flag = 1;
else
    test_flag = 0;
end

% 制图
if plot_flag == 1
    figure(1)
    set(figure(1), 'Position', [50, 200, 1280, 480])
    [x, y] = meshgrid(rho_mesh, z_mesh);
end

% 生成数据协方差矩阵
% len_rho = length(rhoa_obs_log);
% len_phs = length(phs_obs);
% Cd_rho_part = generate_cov(len_rho, log10(f_rhoa), f_Cd);
% Cd_phs_part = generate_cov(len_phs, log10(f_phs), f_Cd);
% Cd = [Cd_rho_part, zeros(len_rho, len_phs); zeros(len_phs, len_rho), Cd_phs_part];
Cd = diag([rhoa_obs_err_log.^2; phs_obs_err.^2]);

% 生成容差系数向量
% k_err_vec_rhoa = k_err .* exp(- rhoa_obs_err_log.^2 ./ (2*((k_err_rhoa.*rhoa_obs_log)/3).^2));
% k_err_vec_phs = k_err .* exp(- phs_obs_err.^2 / (2*(k_err_phs/3)^2));
% k_err_vec = [k_err_vec_rhoa; k_err_vec_phs];

% 处理网格
z_mesh_log = log10(z_mesh);
rho_mesh_log = log10(rho_mesh);

mesh_layers_n = length(z_mesh_log); % 网格层数
mesh_rho_n = length(rho_mesh_log); % 电阻网格数
model_grid = zeros(mesh_layers_n, mesh_rho_n); % 后验概率密度矩阵

% BOSITICK反演生成初始模型
[m_bostick, z_bostick] = bostick_func(rhoa_obs_log, interp1(log10(f_phs), phs_obs, log10(f_rhoa)), f_rhoa); % bostick反演
n = 10; % 初始层数
n_range = [3, 30]; % 层数范围
z_log = [mesh_func(linspace(2, 4, n-1), z_mesh_log); inf];
rho_log = [interp1(z_bostick, m_bostick, z_log(1:end-1)); 0];
rho_log(end) = rho_log(end-1);
rho_log = mesh_func(rho_log, rho_mesh_log);

% 初始化参数
end_flag = 0;
cnt_rejections = 0;
operation = 0;
model_ind = 1;
model_average_log_std = NaN;

% lh = likelihood_func(Cd, rhoa_obs_log, phs_obs, f_rhoa, f_phs, rho_log, z_log, k_weight, k_err_vec); % 计算初始模型的似然函数
lh = likelihood_func(Cd, rhoa_obs_log, phs_obs, f_rhoa, f_phs, rho_log, z_log, k_weight, scale_factor);

Cm = generate_cov(n, z_log, z_smooth_log); % 模型协方差
prior_probability = exp(-norm((Cm^-(1/2))*rho_log)^2 * lambda/2); % 先验（仅包含平滑）
% prior_probability = exp(-norm((Cm^-(1/2))*rho_log)^2/2);
% prior_probability = 1;
ppd = lh * prior_probability;

model_cell = cell(N, 4);
model_cell{1, 1} = z_log; % 分层向量
model_cell{1, 2} = rho_log; % 电阻率向量
model_cell{1, 3} = ppd; % PPD
model_cell{1, 4} = cnt_rejections; % 拒绝次数
model_cell{1, 5} = n; % 层数
model_cell{1, 6} = operation; % 操作类型
model_average_log_matrix = zeros(ceil(N/N_refresh), length(z_mesh));
refresh_ind = 0;

t_func = tic;
while model_ind < N
    flag = rand;
    rho_log_new = rho_log;
    z_log_new = z_log;
    n_new = n;
    
    if flag < 0.34 % 扰动电阻率
%     if mod(model_ind, 2) == 0
        n_rnd = randi(n_new);
        rho_log_new(n_rnd) = perturb_func(rho_log_new(n_rnd), 'rho', rho_mesh_log);
        operation = 1;
        if min(rho_log_new) < rho_mesh_log(2) || max(rho_log_new) > rho_mesh_log(end-1)
            continue;
        end
        lh = likelihood_func(Cd, rhoa_obs_log, phs_obs, f_rhoa, f_phs, rho_log_new, z_log_new, k_weight, scale_factor);
    elseif flag < 0.67 % 扰动分层位置
%     elseif flag < 0.34
        n_rnd = randi(n_new-1);
        z_log_new(n_rnd) = perturb_func(z_log_new(n_rnd), 'z', z_mesh_log);
        operation = 2;
        if z_log_new(1) < z_mesh_log(2) || min(z_log_new(2:end)-z_log_new(1:end-1)) < z_mesh_log(2)
            continue;
        end
        lh = likelihood_func(Cd, rhoa_obs_log, phs_obs, f_rhoa, f_phs, rho_log_new, z_log_new, k_weight, scale_factor);
    else % 生灭
        n_rnd = randi(n_new-1);
        if (rand < 0.5 || n_new >= n_range(2)) && n_new > n_range(1) % 灭
            z_log_new(n_rnd) = [];
            rho_log_new(n_rnd) = [];
            operation = 3;
            n_new = n_new - 1;
            lh = likelihood_func(Cd, rhoa_obs_log, phs_obs, f_rhoa, f_phs, rho_log_new, z_log_new, k_weight, scale_factor) * (2-k_punish);
        else % 生，使第一层和最后一层可以生
            if randi(n_new) == n_new % 使最后一层可以生
                z_log_new = [z_log_new(1:end-1); mesh_func(mean([z_log_new(end-1), z_mesh_log(end)]), z_mesh_log); z_log_new(end)];
%                 z_log_new = [z_log_new(1:end-1); mesh_func(rand * (z_mesh_log(end)-z_log_new(end-1)), z_mesh_log); z_log_new(end)];
                rho_log_new = [rho_log_new(1:end-1); mesh_func((rand * (rho_mesh_log(end)-rho_mesh_log(1))+rho_mesh_log(1)), rho_mesh_log); rho_log_new(end)];
            elseif n_rnd > 1
                z_log_new = [z_log_new(1:n_rnd-1); mesh_func(mean(z_log_new(n_rnd-1:n_rnd)), z_mesh_log); z_log_new(n_rnd:end)];
%                 z_log_new = [z_log_new(1:n_rnd-1); mesh_func(rand * (z_log_new(n_rnd)-z_log_new(n_rnd-1)), z_mesh_log); z_log_new(n_rnd:end)];
                rho_log_new = [rho_log_new(1:n_rnd-1); mesh_func((rand * (rho_mesh_log(end)-rho_mesh_log(1))+rho_mesh_log(1)), rho_mesh_log); rho_log_new(n_rnd:end)];
            else % 使第一层可以生
                z_log_new = [mesh_func(mean([0, z_log_new(n_rnd)]), z_mesh_log); z_log_new(n_rnd:end)];
%                 z_log_new = [mesh_func(rand * z_mesh_log(n_rnd), z_mesh_log); z_log_new(n_rnd:end)];
                rho_log_new = [mesh_func((rand * (rho_mesh_log(end)-rho_mesh_log(1))+rho_mesh_log(1)), rho_mesh_log); rho_log_new(n_rnd:end)];
            end
            operation = 4;
            n_new = n_new + 1;
        lh = likelihood_func(Cd, rhoa_obs_log, phs_obs, f_rhoa, f_phs, rho_log_new, z_log_new, k_weight, scale_factor) * k_punish;
        end
    end
    
    Cm_new = generate_cov(n_new, z_log_new, z_smooth_log);
    prior_probability = exp(-norm((Cm_new^(-1/2))*rho_log_new)^2 * lambda/2);
%     prior_probability_new = 1;
    ppd_new = lh * prior_probability;
    
    if rand < ppd_new/ppd && length(z_log_new) == length(unique(z_log_new)) && isequal(z_log_new, sort(z_log_new)) && max(rho_log_new) < rho_mesh(end)
        model_ind = model_ind + 1;
        model_cell{model_ind, 1} = z_log_new;
        model_cell{model_ind, 2} = rho_log_new;
        model_cell{model_ind, 3} = ppd_new;
        model_cell{model_ind, 4} = cnt_rejections;
        model_cell{model_ind, 5} = n_new;
        model_cell{model_ind, 6} = operation;
        rho_log = rho_log_new;
        z_log = z_log_new;
%         lh = lh_new;
        ppd = ppd_new;
        n = n_new;
        cnt_rejections = 0;
        
        for mesh_n_ind = 1:mesh_layers_n % 计算后验概率
            zk = z_log - z_mesh_log(mesh_n_ind);
            z_ind = find(zk >= 0, 1);
            m_ind = find(rho_mesh_log == rho_log(z_ind));
            model_grid(mesh_n_ind, m_ind) = model_grid(mesh_n_ind, m_ind) + ppd;
        end
        
        if mod(model_ind, N_refresh) == 0
            refresh_ind = refresh_ind + 1;
            % 计算峰值和期望
            rho_average_log = sum(repmat(rho_mesh_log', length(z_mesh_log), 1) .* (model_grid ./ repmat(sum(model_grid, 2), 1, length(rho_mesh_log))), 2);
            [rhoa_average_log, ~] = forward_func(rho_average_log, z_mesh_log, f_rhoa);
            [~, phs_average] = forward_func(rho_average_log, z_mesh_log, f_phs);
            
            [~, rho_max_log_ind] = max(model_grid, [], 2);
            rho_max_log = rho_mesh_log(rho_max_log_ind);
            [rhoa_max_log, ~] = forward_func(rho_max_log, z_mesh_log, f_rhoa);
            [~, phs_max] = forward_func(rho_max_log, z_mesh_log, f_phs);
            
            % 判定期望的变化幅度
            model_average_log_matrix(refresh_ind, :) = rho_average_log';
            if model_ind >= N_end
                model_average_log_std = max(std(model_average_log_matrix(refresh_ind:-1:refresh_ind-round(N_end/N_refresh)+1, :)));
            end
            
            rho_average = 10.^rho_average_log;
            rho_max = 10.^rho_max_log;
            rms_max = rms(([rhoa_obs_log; phs_obs] - [rhoa_max_log; phs_max])./[rhoa_obs_err_log; phs_obs_err]);
            rms_average = rms(([rhoa_obs_log; phs_obs] - [rhoa_average_log; phs_average])./[rhoa_obs_err_log; phs_obs_err]);
            
            % 制图
            if plot_flag == 1
                clf(figure(1))
                annotation('textbox',[.9 0 .1 .9], 'String',['迭代数：', num2str(model_ind)],'EdgeColor','none');
                annotation('textbox',[.9 0 .1 .8], 'String','峰值RMS误差：','EdgeColor','none');
                annotation('textbox',[.9 0 .1 .7], 'String',num2str(rms_max),'EdgeColor','none');
                annotation('textbox',[.9 0 .1 .6], 'String','期望RMS误差：','EdgeColor','none');
                annotation('textbox',[.9 0 .1 .5], 'String',num2str(rms_average),'EdgeColor','none');
                annotation('textbox',[.9 0 .1 .4], 'String',['过去', num2str(N_end), '次迭代'],'EdgeColor','none');
                annotation('textbox',[.9 0 .1 .3], 'String','标准差最大值：','EdgeColor','none');
                annotation('textbox',[.9 0 .1 .2], 'String',num2str(model_average_log_std),'EdgeColor','none');
                
                subplot(1, 4, 2); % 模型
                stairs(rho_average, z_mesh, 'LineWidth', 1)
                set(gca, 'XScale', 'log');
                set(gca, 'YScale', 'log');
                set(gca, 'YDir', 'reverse');
                axis([rho_mesh(1), rho_mesh(end), z_mesh(1), z_mesh(end)]);
                grid on
                hold on
                stairs(rho_max, z_mesh, 'r', 'LineWidth', 1)
                if test_flag == 1
                    plot([m_test(1), m_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
                    stairs(m_test, z_test, 'k--', 'LineWidth', 1);
                    plot([m_test(end), m_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
                    
                    legend('PPD期望', 'PPD峰值', '测试模型')
                else
                    legend('PPD期望', 'PPD峰值')
                end
                xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
                title('PPD峰值与期望')
                xlabel('电阻率')
                ylabel('深度')
                hold off
                
                subplot(1, 4, 1); % 伪彩色图
                mesh(x, y, model_grid, 'FaceColor', 'flat')
                set(gca, 'XScale', 'log');
                set(gca, 'YScale', 'log');
                set(gca, 'ColorScale', 'log');
                view(0,270)
                xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
                title('PPD伪彩色图')
                xlabel('电阻率')
                
                subplot(1, 4, 3); % 正演视电阻率
                semilogy(rhoa_average_log, f_rhoa, 'LineWidth', 1)
                axis([log10(rho_mesh(1)), log10(rho_mesh(end)), f_rhoa(1), f_rhoa(end)]);
                grid on
                hold on
                semilogy(rhoa_max_log, f_rhoa, 'r', 'LineWidth', 1)
                errorbar(rhoa_obs_log, f_rhoa, rhoa_obs_err_log, 'horizontal', 'ko')
                xticks(linspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
                xticklabels(['{10^0}'; '{10^1}'; '{10^2}'; '{10^3}'; '{10^4}'; '{10^5}'; '{10^6}'])
                legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
                title('PPD峰值与期望正演响应')
                xlabel('视电阻率')
                ylabel('频率')
                hold off
                
                subplot(1, 4, 4); % 正演相位
                semilogy(phs_average, f_phs, 'LineWidth', 1)
                axis([0, 90, f_phs(1), f_phs(end)]);
                grid on
                hold on
                semilogy(phs_max, f_phs, 'r', 'LineWidth', 1)
                errorbar(phs_obs, f_phs, phs_obs_err, 'horizontal', 'ko')
                legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
                title('PPD峰值与期望正演响应')
                xlabel('相位')
                ylabel('频率')
                hold off
                
                drawnow;
                % 输出帧
                if print_flag == 1
                    print(gcf, '-djpeg', ['.\Frames\tmp\', num2str(round(model_ind/N_refresh)), '.jpg']);
                end
            end
            
            disp(['迭代数：', num2str(model_ind), '，PPD峰值RMS误差：', num2str(rms_max), '，PPD期望RMS误差：', num2str(rms_average), '，过去', num2str(N_end), '次迭代的模型标准差最大值：', num2str(model_average_log_std)]);
            
            % 判断模型是否已稳定
            if model_average_log_std < std_target && model_ind > N_end
                % 判断模型是否已达标
                t = toc(t_func);
                if (rms_max < rms_target && rms_average < rms_target) || model_ind >= N
                    disp(['迭代结束，本次迭代', num2str(model_ind), '次，耗时', num2str(t), 's']);
                    break;
                end
                disp(['模型已稳定，但不达标，迭代终止。本次迭代', num2str(model_ind), '次，耗时', num2str(t), 's']);
                end_flag = 1;
                break;
            end
        end
        
    else
        cnt_rejections = cnt_rejections + 1;
    end
end