% z_range = [z_mesh(1), z_mesh(end)]; % 坐标轴范围
z_range = [1E1, 1E4];
% rho_range = [rho_mesh(1), rho_mesh(end)];
rho_range = [1E0, 1E4];
test_flag = 0;

%% 分布四联图
figure(5)
set(figure(5), 'Position', [50, 200, 860, 240])

subplot(1, 4, 1); % 伪彩色图
[x, y] = meshgrid(rho_mesh, z_mesh);
mesh(x, y, model_grid_tmp, 'FaceColor', 'flat')
ylim(z_range)
xlim(rho_range)
hold on
grid on
if test_flag == 1
    plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
    stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
    plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
end
stairs(model_conf_edge(:, 3), z_mesh, 'r--', 'LineWidth', 1)
stairs(model_conf_edge(:, 4), z_mesh, 'r--', 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
% set(gca, 'ColorScale', 'log'); % 对数色标
c = colorbar('position', [0.06 0.1 0.02 0.83], 'FontSize', 10);
c.Label.String = 'Probability';
c.Label.FontSize = 12;
view(0,270)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% title('PPD伪彩色图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(1, 4, 2); % 模型
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'YDir', 'reverse');
ylim(z_range)
xlim(rho_range)
grid on
hold on
if test_flag == 1
    plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
    stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
    plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
end
% stairs(rho_median, z_mesh, 'g', 'LineWidth', 1)
stairs(rho_max, z_mesh, 'r', 'LineWidth', 1)
stairs(rho_average, z_mesh, 'b', 'LineWidth', 1)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(1, 4, 3); % 层面概率图
hold on
for j = 2:z_n
    fill([0, z_vec(j), z_vec(j), 0], [z_mesh(j), z_mesh(j), z_mesh(j-1), z_mesh(j-1)], [0.8 0.8 0.8], 'EdgeColor', 'none')
end
if test_flag == 1
    for j = 1:length(z_test)-1
        yline(z_test(j), 'k--', 'LineWidth', 1)
    end
end
stairs(z_vec, z_mesh, 'k');
ylim(z_range)
set(gca, 'YDir', 'reverse');
set(gca, 'YScale', 'log');
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Probability')
ylabel('Interface Depth/m')
grid on


% figure(2)
% set(figure(2), 'Position', [50, 200, 640, 240])
subplot(1, 4, 4); % 层数概率图
[~, peaks_loc] = findpeaks(model_n_hist(:, 2));
hold on
for j = 1:length(model_n_hist)-1
    fill([model_n_hist(j, 1), model_n_hist(j, 1), model_n_hist(j+1, 1),model_n_hist(j+1, 1)], [0, model_n_hist(j, 2), model_n_hist(j, 2), 0], [0.8 0.8 0.8], 'EdgeColor', 'none')
end
stairs(model_n_hist(:, 1), model_n_hist(:, 2), 'k');
% xline(model_n_hist(peaks_loc(1), 1), 'k--', 'LineWidth', 1)
text(model_n_hist(peaks_loc(1), 1),model_n_hist(peaks_loc(1), 2)+0.01,num2str(model_n_hist(peaks_loc(1), 1)),'Color','k', 'FontName','Times New Roman')
grid on
axis([2, 30, 0, 0.2])
% title('模型维数概率密度柱状图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Number of interface')
ylabel('Probability')
xticks(5:5:25)

%% 分布四联图（MCMC）
figure(5)
set(figure(5), 'Position', [50, 200, 860, 240])

subplot(1, 4, 1); % 伪彩色图
[x, y] = meshgrid(rho_mesh, z_mesh);
mesh(x, y, model_grid_tmp, 'FaceColor', 'flat')
ylim(z_range)
xlim(rho_range)
hold on
grid on
if test_flag == 1
    plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
    stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
    plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
end
stairs(model_conf_edge(:, 3), z_mesh, 'r--', 'LineWidth', 1)
stairs(model_conf_edge(:, 4), z_mesh, 'r--', 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
% set(gca, 'ColorScale', 'log'); % 对数色标
c = colorbar('position', [0.06 0.1 0.02 0.83], 'FontSize', 10);
c.Label.String = 'Probability';
c.Label.FontSize = 12;
view(0,270)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% title('PPD伪彩色图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(1, 4, 2); % 模型
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'YDir', 'reverse');
ylim(z_range)
xlim(rho_range)
grid on
hold on
if test_flag == 1
    plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
    stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
    plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
end
% stairs(rho_median, z_mesh, 'g', 'LineWidth', 1)
stairs(rho_max, z_mesh, 'r', 'LineWidth', 1)
stairs(rho_average, z_mesh, 'b', 'LineWidth', 1)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(1, 4, 3); % 视电阻率
semilogy(rhoa_average_log, f_rhoa, 'b', 'LineWidth', 1)
% ylim([f_rhoa(1), f_rhoa(end)]);
ylim([1E-3, 1E3])
xlim([0, 4])
% xlim(rho_range);
grid on
hold on
semilogy(rhoa_max_log, f_rhoa, 'r', 'LineWidth', 1)
errorbar(rhoa_obs_log, f_rhoa, rhoa_obs_err_log, 'horizontal', 'ko')
xticks(linspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
xticklabels(['{10^0}'; '{10^1}'; '{10^2}'; '{10^3}'; '{10^4}'; '{10^5}'; '{10^6}'])
yticks(logspace(log10(f_obs(1)), log10(f_obs(end)), log10(f_obs(end))*2+1))
% legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
% title('PPD峰值与期望正演响应')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Apparent Resistivity/\Omega·m')
ylabel('Frequency/Hz')
hold off


% figure(2)
% set(figure(2), 'Position', [50, 200, 640, 240])
subplot(1, 4, 4); % 相位
semilogy(phs_average, f_phs, 'b', 'LineWidth', 1)
xlim([0, 90]);
ylim([1E-3, 1E3])
grid on
hold on
semilogy(phs_max, f_phs, 'r', 'LineWidth', 1)
errorbar(phs_obs, f_phs, phs_obs_err, 'horizontal', 'ko')
yticks(logspace(log10(f_obs(1)), log10(f_obs(end)), log10(f_obs(end))*2+1))
% legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
% title('PPD峰值与期望正演响应')
xticks([0, 30, 60, 90])
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Phase/°')
ylabel('Frequency/Hz')
hold off
%% 分布六联图
f_obs = f_rhoa;

figure(1)
% set(figure(1), 'Position', [0, 0, 640, 480])
set(figure(1), 'Position', [0, 0, 640, 720])

subplot(2, 3, 1); % 伪彩色图
[x, y] = meshgrid(rho_mesh, z_mesh);
mesh(x, y, model_grid_tmp, 'FaceColor', 'flat')
ylim(z_range)
% xlim([1E0, 1E4])
xlim(rho_range);
hold on
grid on
if test_flag == 1
    plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
    stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
    plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
end
stairs(model_conf_edge(:, 3), z_mesh, 'r--', 'LineWidth', 1)
stairs(model_conf_edge(:, 4), z_mesh, 'r--', 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
% set(gca, 'ColorScale', 'log'); % 对数色标
c = colorbar('position', [0.05 0.59 0.02 0.33], 'FontSize', 10);
c.Label.String = 'Probability';
c.Label.FontSize = 12;
view(0,270)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
yticks(logspace(log10(z_mesh(1)), log10(z_mesh(end)), log10(z_mesh(end))+1))
% title('PPD伪彩色图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(2, 3, 2); % 模型
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'YDir', 'reverse');
ylim(z_range)
% xlim([1E0, 1E4])
xlim(rho_range);
grid on
hold on
if test_flag == 1
plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
end
stairs(rho_average, z_mesh, 'b', 'LineWidth', 1)
% stairs(rho_median, z_mesh, 'g', 'LineWidth', 1)
stairs(rho_max, z_mesh, 'r', 'LineWidth', 1)
stairs(rho_occam, z_mesh_occam, 'g', 'LineWidth', 1) % OCCAM
% plot(occam_test(:, 2), occam_test(:, 1), 'm', 'LineWidth', 1) % OCCAM_test

xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
yticks(logspace(log10(z_mesh(1)), log10(z_mesh(end)), log10(z_mesh(end))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(2, 3, 3); % 层面概率图
hold on
for j = 2:z_n
    fill([0, z_vec(j), z_vec(j), 0], [z_mesh(j), z_mesh(j), z_mesh(j-1), z_mesh(j-1)], [0.8 0.8 0.8], 'EdgeColor', 'none')
end
% for j = 1:length(z_test)-1
%     yline(z_test(j), 'k--', 'LineWidth', 1)
% end
stairs(z_vec, z_mesh, 'k');
ylim(z_range)
set(gca, 'YDir', 'reverse');
set(gca, 'YScale', 'log');
set(gca,'FontName','Times New Roman', 'FontSize', 10)
yticks(logspace(log10(z_mesh(1)), log10(z_mesh(end)), log10(z_mesh(end))+1))
xlabel('Probability')
ylabel('Interface Depth/m')
grid on


% figure(2)
% set(figure(2), 'Position', [50, 200, 640, 240])
subplot(2, 3, 4); % 层数概率图
[~, peaks_loc] = findpeaks(model_n_hist(:, 2));
hold on
for j = 1:length(model_n_hist)-1
    fill([model_n_hist(j, 1), model_n_hist(j, 1), model_n_hist(j+1, 1),model_n_hist(j+1, 1)], [0, model_n_hist(j, 2), model_n_hist(j, 2), 0], [0.8 0.8 0.8], 'EdgeColor', 'none')
end
stairs(model_n_hist(:, 1), model_n_hist(:, 2), 'k');
% xline(model_n_hist(peaks_loc(1), 1), 'k--', 'LineWidth', 1)
text(model_n_hist(peaks_loc(1), 1),model_n_hist(peaks_loc(1), 2)+0.01,num2str(model_n_hist(peaks_loc(1), 1)),'Color','k',  'FontName','Times New Roman')
grid on
axis([2, 30, 0, 0.2])
% title('模型维数概率密度柱状图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Number of interface')
ylabel('Probability')
xticks(5:5:25)

subplot(2, 3, 5); % 正演视电阻率
semilogy(rhoa_average_log, f_rhoa, 'b', 'LineWidth', 1)
% ylim([f_rhoa(1), f_rhoa(end)]);
ylim([10^floor(log10(f_obs(1))), 10^ceil(log10(f_obs(end)))])
xlim([-1, 3])
% xlim(rho_range);
grid on
hold on
semilogy(rhoa_max_log, f_rhoa, 'r', 'LineWidth', 1)
semilogy(rhoa_occam_log, f_rhoa, 'g', 'LineWidth', 1)
% semilogy(rhoa_occam_test(:, 2), rhoa_occam_test(:, 1), 'm', 'LineWidth', 1)
errorbar(rhoa_obs_log, f_rhoa, rhoa_obs_err_log, 'horizontal', 'ko')
% xticks(linspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
xticklabels({'{10^{-1}}'; '{10^0}'; '{10^1}'; '{10^2}'; '{10^3}'; '{10^4}'; '{10^5}'; '{10^6}'})
yticks(logspace(floor(log10(f_obs(1))), ceil(log10(f_obs(end))), ceil(log10(f_obs(end)))-floor(log10(f_obs(1)))+1))
% legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
% title('PPD峰值与期望正演响应')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Apparent Resistivity/\Omega·m')
ylabel('Frequency/Hz')
hold off

subplot(2, 3, 6); % 正演相位
semilogy(phs_average, f_phs, 'b', 'LineWidth', 1)
xlim([0, 90]);
ylim([10^floor(log10(f_obs(1))), 10^ceil(log10(f_obs(end)))])
grid on
hold on
semilogy(phs_max, f_phs, 'r', 'LineWidth', 1)
errorbar(phs_obs, f_phs, phs_obs_err, 'horizontal', 'ko')
semilogy(phs_occam, f_obs, 'g', 'LineWidth', 1)
yticks(logspace(floor(log10(f_obs(1))), ceil(log10(f_obs(end))), ceil(log10(f_obs(end)))-floor(log10(f_obs(1)))+1))
% legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
% title('PPD峰值与期望正演响应')
xticks([0, 30, 60, 90])
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Phase/°')
ylabel('Frequency/Hz')
hold off

%% 误差
figure(2)
set(figure(2), 'Position', [50, 200, 640, 240])
loglog([model_cell{:, 7}])
hold on
grid on
xline(N_iter_burn_in, 'k--', 'LineWidth', 1)
text(N_iter_burn_in-7E3,1E1,'Burn-in','Color','k', 'FontSize', 12,'FontName','Times New Roman')
% axis([1, 1E5, 3E1, 2E3]);
set(gca,'FontName','Times New Roman', 'FontSize', 10)
ylabel('Data Square Error')
xlabel('Sample')

%% 图例
figure(3)
ax1 = axes('Position',[0.13 0.58 0, 0]);
yline(1, 'r--', 'LineWidth', 1)
hold on
yline(1, 'k--', 'LineWidth', 1)
yline(1, 'r', 'LineWidth', 1)
yline(1, 'b', 'LineWidth', 1)
% yline(1, 'g', 'LineWidth', 1)
errorbar(0, 1, 1, 'horizontal', 'ko')
% legend('95% Confidence Interval', 'True Model', 'Maximum Probability Model', 'Mathematical Expected Model', 'Observation Data', 'NumColumns',3, 'EdgeColor', 'none')
legend('95% Confidence Interval', 'True Model', 'Maximum Probability Model', 'Mathematical Expected Model', 'NumColumns',3, 'EdgeColor', 'none')
set(gca,'FontName','Times New Roman', 'FontSize', 10)

figure(8)
ax1 = axes('Position',[0.13 0.58 0, 0]);
yline(1, 'r--', 'LineWidth', 1)
hold on
yline(1, 'k--', 'LineWidth', 1)
% plot(1, 1, 'kx')
yline(1, 'r', 'LineWidth', 1)
yline(1, 'b', 'LineWidth', 1)
% yline(1, 'g', 'LineWidth', 1)
errorbar(0, 1, 1, 'horizontal', 'ko')
legend('95% Confidence Interval', 'True Model', 'Maximum Probability Model', 'Mathematical Expected Model', 'Observation Data', 'NumColumns',3, 'EdgeColor', 'none')
% legend('95% Confidence Interval', 'Normal Resistivity Log Data', 'Maximum Probability Model', 'Mathematical Expected Model', 'OCCAM Model', 'Observation Data', 'NumColumns',3, 'EdgeColor', 'none')
set(gca,'FontName','Times New Roman', 'FontSize', 10)

figure(9)
ax1 = axes('Position',[0.13 0.58 0, 0]);
yline(1, 'r--', 'LineWidth', 1)
hold on
% yline(1, 'k--', 'LineWidth', 1)
yline(1, 'r', 'LineWidth', 1)
yline(1, 'b', 'LineWidth', 1)
yline(1, 'g', 'LineWidth', 1)
yline(1, 'm', 'LineWidth', 1)
errorbar(0, 1, 1, 'horizontal', 'ko')
% legend('95% Confidence Interval', 'True Model', 'Maximum Probability Model', 'Mathematical Expected Model', 'Observation Data', 'NumColumns',3, 'EdgeColor', 'none')
legend('95% Confidence Interval', 'Maximum Probability Model', 'Mathematical Expected Model', 'OCCAM Model', 'M.T. OCCAM Model', 'Observation Data', 'NumColumns',3, 'EdgeColor', 'none')
set(gca,'FontName','Times New Roman', 'FontSize', 10)

%% OCCAM
figure(4)
set(figure(4), 'Position', [50, 200, 640, 240])
z_mesh_occam = logspace(0,5,50)';
z_mesh_occam_mod = [z_mesh_occam(1, 1); z_mesh_occam(2:end)-z_mesh_occam(1:end-1)];
z_mesh_occam_mod = [z_mesh_occam_mod, 2*ones(50, 1)];

subplot(1, 3, 1)
stairs(rho_occam, z_mesh_occam, 'g', 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'YDir', 'reverse');
ylim(z_range)
xlim(rho_range)
grid on
hold on
% plot(occam_test(:, 2), occam_test(:, 1), 'm', 'LineWidth', 1) % OCCAM_test
if test_flag == 1
    plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
    stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
    plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
end
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% legend('光滑反演', '测试模型')
% title('光滑反演')
yticks(logspace(log10(z_mesh(1)), log10(z_mesh(end)), log10(z_mesh(end))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')
% plot(1000./point_name(:, 2), point_name(:, 1), 'kx', 'LineWidth', 1)

[rhoa_occam_log, phs_occam] = forward_func(log10(rho_occam(1:end-1))', log10(z_mesh_occam), f_obs);
subplot(1, 3, 2)
semilogy(rhoa_occam_log, f_obs, 'g', 'LineWidth', 1)
ylim([10^floor(log10(f_obs(1))), 10^ceil(log10(f_obs(end)))])
% xlim(log10(rho_range))
xlim([0, 4])
hold on
grid on
% semilogy(rhoa_occam_test(:, 2), rhoa_occam_test(:, 1), 'm', 'LineWidth', 1)
errorbar(rhoa_obs_log, f_obs, rhoa_obs_err_log, 'horizontal', 'ko')
xticks(linspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
xticklabels({'{10^0}'; '{10^1}'; '{10^2}'; '{10^3}'; '{10^4}'; '{10^5}'; '{10^6}'})
yticks(logspace(floor(log10(f_obs(1))), ceil(log10(f_obs(end))), ceil(log10(f_obs(end)))-floor(log10(f_obs(1)))+1))
% legend('光滑反演正演响应', '观测数据')
% title('光滑反演正演响应')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Apparent Resistivity/\Omega·m')
ylabel('Frequency/Hz')

subplot(1, 3, 3); % 正演相位
semilogy(phs_occam, f_obs, 'g', 'LineWidth', 1)
xlim([0, 90])
ylim([10^floor(log10(f_obs(1))), 10^ceil(log10(f_obs(end)))])
grid on
hold on
% semilogy(phs_occam_test(:, 2), phs_occam_test(:, 1), 'm', 'LineWidth', 1)
errorbar(phs_obs, f_obs, phs_obs_err, 'horizontal', 'ko')
% legend('光滑反演正演响应', '观测数据')
% title('光滑反演正演响应')
xticks([0, 30, 60, 90])
yticks(logspace(floor(log10(f_obs(1))), ceil(log10(f_obs(end))), ceil(log10(f_obs(end)))-floor(log10(f_obs(1)))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Phase/°')
ylabel('Frequency/Hz')
hold off

%% 实测数据制图
% MK6
H08_res = [29.7, 570; 42.5, 70; 45, 70; 134.5, 70; 603.4, 36; 960.75, 110; 970, 35; 985.2, 60; 1026.7, 52; 1332.5, 82];

% MK2
F06_res = [48.9, 50; 301, 10; 673.7, 120; 768.1, 10; 833, 35; 895.5, 55; 960.2, 70; 1082.3, 72; 1084.6, 85];
point_name = F06_res;
rho_range = [1E0, 1E4];
z_range = [1E0, 1E4];

figure(6)
set(figure(6), 'Position', [0, 0, 640, 720])

subplot(2, 3, 1); % 伪彩色图
[x, y] = meshgrid(rho_mesh, z_mesh);
mesh(x, y, model_grid_tmp, 'FaceColor', 'flat')
ylim(z_range)
% xlim([1E0, 1E4])
xlim(rho_range);
hold on
grid on
stairs(model_conf_edge(:, 3), [z_mesh(2:end); 10.^(z_mesh_log(end)+z_mesh_log(2)-z_mesh_log(1))], 'r--', 'LineWidth', 1)
stairs(10.^(log10(model_conf_edge(:, 4))+rho_mesh_log(2)-rho_mesh_log(1)), z_mesh, 'r--', 'LineWidth', 1)
plot(1000./point_name(:, 2), point_name(:, 1), 'kx', 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
% set(gca, 'ColorScale', 'log'); % 对数色标
c = colorbar('position', [0.05 0.59 0.02 0.33], 'FontSize', 10);
c.Label.String = 'Probability';
c.Label.FontSize = 12;
view(0,270)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% title('PPD伪彩色图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(2, 3, 2); % 模型
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'YDir', 'reverse');
ylim(z_range)
% xlim([1E0, 1E4])
xlim(rho_range);
grid on
hold on
stairs(rho_occam, z_mesh_occam, 'g', 'LineWidth', 1)
stairs(rho_average, z_mesh, 'b', 'LineWidth', 1)
stairs(rho_max, z_mesh, 'r', 'LineWidth', 1)
plot(1000./point_name(:, 2), point_name(:, 1), 'kx', 'LineWidth', 1)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(2, 3, 3); % 层面概率图
hold on
for j = 2:z_n
    fill([0, z_vec(j), z_vec(j), 0], [z_mesh(j), z_mesh(j), z_mesh(j-1), z_mesh(j-1)], [0.8 0.8 0.8], 'EdgeColor', 'none')
end
% for j = 1:length(z_test)-1
%     yline(z_test(j), 'k--', 'LineWidth', 1)
% end
stairs(z_vec, z_mesh, 'k');
ylim(z_range)
set(gca, 'YDir', 'reverse');
set(gca, 'YScale', 'log');
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Probability')
ylabel('Interface Depth/m')
grid on


% figure(2)
% set(figure(2), 'Position', [50, 200, 640, 240])
subplot(2, 3, 4); % 层数概率图
[~, peaks_loc] = findpeaks(model_n_hist(:, 2));
hold on
for j = 1:length(model_n_hist)-1
    fill([model_n_hist(j, 1), model_n_hist(j, 1), model_n_hist(j+1, 1),model_n_hist(j+1, 1)], [0, model_n_hist(j, 2), model_n_hist(j, 2), 0], [0.8 0.8 0.8], 'EdgeColor', 'none')
end
stairs(model_n_hist(:, 1), model_n_hist(:, 2), 'k');
% xline(model_n_hist(peaks_loc(1), 1), 'k--', 'LineWidth', 1)
text(model_n_hist(peaks_loc(1), 1),model_n_hist(peaks_loc(1), 2)+0.01,num2str(model_n_hist(peaks_loc(1), 1)),'Color','k',  'FontName','Times New Roman')
grid on
axis([2, 30, 0, 0.2])
% title('模型维数概率密度柱状图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Number of interface')
ylabel('Probability')
xticks(5:5:25)

subplot(2, 3, 5); % 正演视电阻率
ylim([f_rhoa(1), f_rhoa(end)]);
% ylim([1E-3, 1E3])
xlim([0, 4])
% xlim(rho_range);
grid on
hold on
semilogy(rhoa_occam_log, f_obs, 'g', 'LineWidth', 1)
semilogy(rhoa_average_log, f_rhoa, 'b', 'LineWidth', 1)
semilogy(rhoa_max_log, f_rhoa, 'r', 'LineWidth', 1)
errorbar(rhoa_obs_log, f_rhoa, rhoa_obs_err_log, 'horizontal', 'ko')
xticks(linspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
xticklabels(['{10^0}'; '{10^1}'; '{10^2}'; '{10^3}'; '{10^4}'; '{10^5}'; '{10^6}'])
% yticks(logspace(log10(f_obs(1)), log10(f_obs(end)), log10(f_obs(end))*2+1))
% legend('PPD期望正演响应', 'PPD峰值正演响应', '观测数据')
% title('PPD峰值与期望正演响应')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
set(gca, 'YScale', 'log');
xlabel('Apparent Resistivity/\Omega·m')
ylabel('Frequency/Hz')
hold off

subplot(2, 3, 6); % 正演相位
xlim([0, 90]);
% ylim([1E-3, 1E3])
ylim([f_rhoa(1), f_rhoa(end)]);
grid on
hold on
semilogy(phs_occam, f_obs, 'g', 'LineWidth', 1)
semilogy(phs_average, f_phs, 'b', 'LineWidth', 1)
semilogy(phs_max, f_phs, 'r', 'LineWidth', 1)
errorbar(phs_obs, f_phs, phs_obs_err, 'horizontal', 'ko')
xticks([0, 30, 60, 90])
set(gca,'FontName','Times New Roman', 'FontSize', 10)
set(gca, 'YScale', 'log');
xlabel('Phase/°')
ylabel('Frequency/Hz')
hold off

%% 岩性柱对比

figure(7)
set(figure(7), 'Position', [0, 50, 320, 720])
% MK6、MK2
I08_Lith = [0, 0; 40, 0; 70, 1; 100, 2; 360, 3; 400, 1; 410, 5; 450, 3; 480, 4; 500, 1; 800, 4; 810, 5; 990, 4; 1000, 5; 1090, 4; 1700, 6];
G05_Lith = [0, 0; 30, 0; 60, 1; 70, 2; 150, 0; 380, 1; 390, 3; 410, 1; 415, 4; 420, 5; 490, 7; 520, 6; 680, 5; 1200, 8];
% point_Lith = I08_Lith;
point_Lith = G05_Lith;

I08_Lith_2 = [0, 0; 70, 0; 450, 1; 480, 2; 500, 1; 1700, 2];
G05_Lith_2 = [0, 0; 70, -1; 150, 0; 380, 1; 390, 2; 420, 1; 490, 2; 680, 1; 1200, 2];
% point_Lith_2 = I08_Lith_2;
point_Lith_2 = G05_Lith_2;

subplot(1, 3, 1)
hold on
set(gca, 'YDir', 'reverse');
set(gca,'tickdir','out')
set(gca,'position', [0.18 0.23 0.05 0.76]);
ylim([point_Lith(1, 1), point_Lith(end, 1)])
c = zeros(9, 1);
for j = 1:length(point_Lith)-1
    switch point_Lith(j+1, 2)
        case 0
            color_name = 'w';
        case 1
            color_name = 'r';
        case 2
            color_name = 'b';
        case 3
            color_name = 'g';
        case 4
            color_name = 'y';
        case 5
            color_name = 'c';
        case 6
            color_name = 'm';
        case 7
            color_name = [0.929, 0.694, 0.125];
        case 8
            color_name = [0.494, 0.184, 0.556];
    end
    hdl(j) = fill([0, 1000, 1000, 0], [point_Lith(j, 1), point_Lith(j, 1), point_Lith(j+1, 1), point_Lith(j+1, 1)], color_name);
    if c(point_Lith(j+1, 2)+1) == 0
        c(point_Lith(j+1, 2)+1) = 1;
    else
        set(hdl(j),'handlevisibility','off');
    end
end
% legend('Non Core', 'Tuff Breccia', 'Andesite', 'Tuff', 'Prophyry', 'Tuff', 'Oz Pyrite', 'NumColumns', 3, 'EdgeColor', 'none', 'Position', [0.4 0.1 0.2 0.05]) % H08
legend('Non Core', 'Dacite', 'Propylite', 'Andesite', 'Breccia', 'Tuff', 'Tuff Breccia', 'Diolite', 'Basalt', 'NumColumns', 3, 'EdgeColor', 'none', 'Position', [0.4 0.1 0.2 0.05]) % F06
set(gca,'FontName','Times New Roman', 'FontSize', 10)
ylabel('Depth/m')
yticks(0:100:point_Lith(end, 1))
xticks([0, 1000])
xticklabels('');

subplot(1, 3, 2)
hold on
set(gca, 'YDir', 'reverse');
% set(gca,'tickdir','out')
set(gca,'position', [0.23 0.23 0.05 0.76]);
ylim([point_Lith(1, 1), point_Lith(end, 1)])
c = zeros(9, 1);
for j = 1:length(point_Lith_2)-1
    hdl(j) = fill([0, 1000, 1000, 0], [point_Lith_2(j, 1), point_Lith_2(j, 1), point_Lith_2(j+1, 1), point_Lith_2(j+1, 1)], 'w');
end
% legend('Non Core', 'Tuff Breccia', 'Andesite', 'Tuff', 'Prophyry', 'Tuff', 'Oz Pyrite', 'NumColumns', 3, 'EdgeColor', 'none', 'Position', [0.4 0.1 0.2 0.05]) % H08
% legend('Non Core', 'Dacite', 'Propylite', 'Andesite', 'Breccia', 'Tuff', 'Tuff Breccia', 'Diolite', 'Basalt', 'NumColumns', 3, 'EdgeColor', 'none', 'Position', [0.4 0.1 0.2 0.05]) % F06
set(gca,'FontName','Times New Roman', 'FontSize', 10)
% ylabel('Depth/m')
yticks(0:100:point_Lith(end, 1))
yticklabels('');
xticks([0, 1000])
xticklabels('');

subplot(1, 3, 3); % 伪彩色图
set(gca,'position', [0.3 0.23 0.6 0.76]);
[x, y] = meshgrid(rho_mesh, z_mesh);
mesh(x, y, model_grid_tmp, 'FaceColor', 'flat', 'handlevisibility','off')
ylim([point_Lith(1, 1), point_Lith(end, 1)])
% xlim([1E0, 1E4])
xlim(rho_range);
hold on
grid on
stairs(model_conf_edge(:, 3), [z_mesh(2:end); 10.^(z_mesh_log(end)+z_mesh_log(2)-z_mesh_log(1))], 'r--', 'LineWidth', 1)
stairs(10.^(log10(model_conf_edge(:, 4))+rho_mesh_log(2)-rho_mesh_log(1)), z_mesh, 'r--', 'LineWidth', 1,'handlevisibility','off');
plot(1000./point_name(:, 2), point_name(:, 1), 'kx', 'LineWidth', 1)
stairs(rho_average, z_mesh, 'w', 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca,'tickdir','out')
% set(gca, 'YScale', 'log');
% set(gca, 'ColorScale', 'log'); % 对数色标
% c = colorbar('position', [0.05 0.59 0.02 0.33], 'FontSize', 10);
% c.Label.String = 'Probability';
% c.Label.FontSize = 12;
view(0,270)
legend('95% Confidence Interval', 'Normal Resistivity Log Data', 'Mathematical Expected Model', 'NumColumns', 1, 'EdgeColor', 'none', 'Position', [0.234 0.04 0.2 0.05]) % F06
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% title('PPD伪彩色图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
yticks(0:100:point_Lith(end, 1))
yticklabels('')
% ylabel('Depth/m')

% 部分标注
% annotation('textbox', [0.66 0.09 0.2 0.05], 'String', {'Ntv:Nakitsurayama', '       Volcanics', 'Si   :Siodomarigawa', '     Formation', 'Int :Intrusion'}, 'FitBoxToText', 'on', 'EdgeColor', 'none','FontName','Times New Roman', 'FontSize', 9);

annotation('textbox', [0.62 0.05 0.2 0.05], 'String', {'Ykv Yokotsu Volcanics', 'Si     Siodomarigawa', '        Formation', 'Int    Intrusion'}, 'FitBoxToText', 'on', 'EdgeColor', 'none','FontName','Times New Roman', 'FontSize', 9);
% annotation('textarrow', [0.5, 0.25], [0.92, 0.92]) % MK2
% annotation('textbox', [0.5 0.888 0.1 0.05], 'String', 'Yokotsu Volcanics', 'FitBoxToText', 'on', 'Color', 'White', 'EdgeColor', 'none','FontName','Times New Roman', 'FontSize', 11);
