

figure(1)
test_flag = 1;

set(figure(1), 'Position', [50, 200, 1280, 480])
[x, y] = meshgrid(rho_mesh, z_mesh);
subplot(1, 4, 1); % 模型
loglog(model_average, z_mesh, 'LineWidth', 1)
set(gca, 'YDir', 'reverse');
axis([rho_mesh(1), rho_mesh(end), z_mesh(1), z_mesh(end)]);
grid on
hold on
loglog(model_max, z_mesh, 'r', 'LineWidth', 1)
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

subplot(1, 4, 2); % 伪彩色图
mesh(x, y, model_grid, 'FaceColor', 'flat')
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'ColorScale', 'log');
view(0,270)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
title('PPD伪彩色图')
xlabel('电阻率')

subplot(1, 4, 3); % 正演视电阻率
semilogy(rho_average_log, f_obs, 'LineWidth', 1)
axis([log10(rho_mesh(1)), log10(rho_mesh(end)), f_obs(1), f_obs(end)]);
grid on
hold on
semilogy(rho_max_log, f_obs, 'r', 'LineWidth', 1)
errorbar(d_obs_log, f_obs, d_obs_err_log, 'horizontal', 'ko')
xticks(linspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
xticklabels(['{10^0}'; '{10^1}'; '{10^2}'; '{10^3}'; '{10^4}'; '{10^5}'; '{10^6}'])
legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
title('PPD峰值与期望正演响应')
xlabel('视电阻率')
ylabel('频率')
hold off

subplot(1, 4, 4); % 正演相位
semilogy(phs_average, f_obs, 'LineWidth', 1)
axis([0, 90, f_obs(1), f_obs(end)]);
grid on
hold on
semilogy(phs_max, f_obs, 'r', 'LineWidth', 1)
errorbar(phs_obs, f_obs, phs_obs_err, 'horizontal', 'ko')
legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
title('PPD峰值与期望正演响应')
xlabel('相位')
ylabel('频率')
hold off

figure(2)
set(figure(2), 'Position', [50, 200, 320, 480])

bar(model_n_hist(:, 2));
grid on
axis([3-model_n_hist(1), 30-model_n_hist(1), 0, 0.35])
title('模型维数概率密度柱状图')
xlabel('模型层数')
ylabel('概率')
xticks(6-model_n_hist(1):5:length(model_n_hist))
xticklabels(num2str(model_n_hist(1:5:end, 1) - (model_n_hist(1, 1) - 5)));

% %% 3
% figure(3)
% set(figure(3), 'Position', [50, 200, 960, 480])
% 
% z_mesh_occam = logspace(-1,5,25)';
% z_mesh_occam_mod = [z_mesh_occam(1, 1); z_mesh_occam(2:end)-z_mesh_occam(1:end-1)];
% z_mesh_occam_mod = [z_mesh_occam_mod, 2*ones(25, 1)];
% 
% subplot(1, 3, 1)
% stairs(10.^resi(1:end-1), z_mesh_occam)
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% set(gca, 'YDir', 'reverse');
% axis([rho_mesh(1), rho_mesh(end), z_mesh(1), z_mesh(end)]);
% grid on
% hold on
% plot([m_test(1), m_test(1)], [z_mesh(1), z_test(1)], 'k--')
% stairs(m_test, z_test, 'k--');
% plot([m_test(end), m_test(end)], [z_test(end-1), z_mesh(end)], 'k--');
% xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% legend('光滑反演', '测试模型')
% title('光滑反演')
% xlabel('电阻率')
% ylabel('深度')
% 
% [occam_rho_log, occam_phs] = forward_func(resi(1:end-1)', log10(z_mesh_occam), f_obs);
% % occam_rho_log = flipud(occam_rho_log');
% % occam_phs = occam_phs';
% subplot(1, 3, 2)
% semilogy(occam_rho_log, f_obs)
% axis([log10(rho_mesh(1)), log10(rho_mesh(end)), f_obs(1), f_obs(end)]);
% hold on
% grid on
% errorbar(d_obs_log, f_obs, d_obs_err_log, 'horizontal', 'ko')
% xticks(linspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% xticklabels(['{10^0}'; '{10^1}'; '{10^2}'; '{10^3}'; '{10^4}'; '{10^5}'; '{10^6}'])
% legend('光滑反演正演响应', '观测数据')
% title('光滑反演正演响应')
% xlabel('视电阻率')
% ylabel('频率')
% 
% subplot(1, 3, 3); % 正演相位
% semilogy(occam_phs, f_obs)
% axis([0, 90, f_obs(1), f_obs(end)]);
% grid on
% hold on
% errorbar(phs_obs, f_obs, phs_obs_err, 'horizontal', 'ko')
% legend('光滑反演正演响应', '观测数据')
% title('光滑反演正演响应')
% xlabel('相位')
% ylabel('频率')
% hold off
