z_range = [1E2, z_mesh(end)]; % 坐标轴范围

%% 分布四联图
figure(5)
set(figure(5), 'Position', [50, 200, 860, 240])

subplot(1, 4, 1); % 伪彩色图
[x, y] = meshgrid(rho_mesh, z_mesh);
mesh(x, y, model_grid_tmp, 'FaceColor', 'flat')
ylim(z_range)
xlim([1E0, 1E4])
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
% c = colorbar('position', [0.05 0.1 0.02 0.83], 'FontSize', 8);
% c.Label.String = 'Probability';
view(0,270)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% title('PPD伪彩色图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(1, 4, 2); % 模型
plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'YDir', 'reverse');
ylim(z_range)
xlim([1E0, 1E4])
grid on
hold on
stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
stairs(rho_average, z_mesh, 'b', 'LineWidth', 1)
% stairs(rho_median, z_mesh, 'g', 'LineWidth', 1)
stairs(rho_max, z_mesh, 'r', 'LineWidth', 1)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(1, 4, 3); % 层面概率图
hold on
for j = 2:z_n
    fill([0, z_vec(j), z_vec(j), 0], [z_mesh(j), z_mesh(j), z_mesh(j-1), z_mesh(j-1)], [0.8 0.8 0.8], 'EdgeColor', 'none')
end
for j = 1:length(z_test)-1
    yline(z_test(j), 'k--', 'LineWidth', 1)
end
stairs(z_vec, z_mesh, 'k');
ylim(z_range)
set(gca, 'YDir', 'reverse');
set(gca, 'YScale', 'log');
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Probability')
ylabel('Depth/m')
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

%% 分布六联图
figure(1)
set(figure(1), 'Position', [0, 0, 640, 480])

subplot(2, 3, 1); % 伪彩色图
[x, y] = meshgrid(rho_mesh, z_mesh);
mesh(x, y, model_grid_tmp, 'FaceColor', 'flat')
ylim(z_range)
xlim([1E0, 1E4])
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
% c = colorbar('position', [0.05 0.1 0.02 0.83], 'FontSize', 8);
% c.Label.String = 'Probability';
view(0,270)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% title('PPD伪彩色图')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(2, 3, 2); % 模型
plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'YDir', 'reverse');
ylim(z_range)
xlim([1E0, 1E4])
grid on
hold on
stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
stairs(rho_average, z_mesh, 'b', 'LineWidth', 1)
% stairs(rho_median, z_mesh, 'g', 'LineWidth', 1)
stairs(rho_max, z_mesh, 'r', 'LineWidth', 1)
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

subplot(2, 3, 3); % 层面概率图
hold on
for j = 2:z_n
    fill([0, z_vec(j), z_vec(j), 0], [z_mesh(j), z_mesh(j), z_mesh(j-1), z_mesh(j-1)], [0.8 0.8 0.8], 'EdgeColor', 'none')
end
for j = 1:length(z_test)-1
    yline(z_test(j), 'k--', 'LineWidth', 1)
end
stairs(z_vec, z_mesh, 'k');
ylim(z_range)
set(gca, 'YDir', 'reverse');
set(gca, 'YScale', 'log');
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Probability')
ylabel('Depth/m')
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
ylim([f_rhoa(1), f_rhoa(end)]);
xlim([0, 4])
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

subplot(2, 3, 6); % 正演相位
semilogy(phs_average, f_phs, 'b', 'LineWidth', 1)
axis([0, 90, f_phs(1), f_phs(end)]);
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

%% 误差
figure(2)
set(figure(2), 'Position', [50, 200, 640, 240])
loglog([model_cell{:, 7}])
hold on
grid on
xline(N_iter_burn_in, 'k--', 'LineWidth', 1)
text(N_iter_burn_in-7E3,1E3,'Burn-in','Color','k')
axis([1, 1E5, 3E1, 2E3]);
set(gca,'FontName','Times New Roman', 'FontSize', 10)
ylabel('Data L2 Loss')
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

%% OCCAM
figure(4)
set(figure(4), 'Position', [50, 200, 640, 240])
z_mesh_occam = logspace(-1,5,100)';
z_mesh_occam_mod = [z_mesh_occam(1, 1); z_mesh_occam(2:end)-z_mesh_occam(1:end-1)];
z_mesh_occam_mod = [z_mesh_occam_mod, 2*ones(100, 1)];

subplot(1, 3, 1)
stairs(10.^resi(1:end-1), z_mesh_occam, 'LineWidth', 1)
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'YDir', 'reverse');
ylim(z_range);
xlim([1E0, 1E4])
grid on
hold on
plot([rho_test(1), rho_test(1)], [z_mesh(1), z_test(1)], 'k--', 'LineWidth', 1)
stairs(rho_test, z_test, 'k--', 'LineWidth', 1);
plot([rho_test(end), rho_test(end)], [z_test(end-1), z_mesh(end)], 'k--', 'LineWidth', 1);
xticks(logspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
% legend('光滑反演', '测试模型')
% title('光滑反演')
yticks(logspace(log10(z_mesh(1)), log10(z_mesh(end)), log10(z_mesh(end))+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Resistivity/\Omega·m')
ylabel('Depth/m')

[occam_rho_log, occam_phs] = forward_func(resi(1:end-1)', log10(z_mesh_occam), f_obs);
% occam_rho_log = flipud(occam_rho_log');
% occam_phs = occam_phs';
subplot(1, 3, 2)
semilogy(occam_rho_log, f_obs, 'LineWidth', 1)
ylim([f_obs(1), f_obs(end)]);
xlim([0, 4])
hold on
grid on
errorbar(rhoa_obs_log, f_obs, rhoa_obs_err_log, 'horizontal', 'ko')
xticks(linspace(log10(rho_mesh(1)), log10(rho_mesh(end)), log10(rho_mesh(end))+1))
xticklabels(['{10^0}'; '{10^1}'; '{10^2}'; '{10^3}'; '{10^4}'; '{10^5}'; '{10^6}'])
yticks(logspace(log10(f_obs(1)), log10(f_obs(end)), log10(f_obs(end))*2+1))
% legend('光滑反演正演响应', '观测数据')
% title('光滑反演正演响应')
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Apparent Resistivity/\Omega·m')
ylabel('Frequency/Hz')

subplot(1, 3, 3); % 正演相位
semilogy(occam_phs, f_obs, 'LineWidth', 1)
axis([0, 90, f_obs(1), f_obs(end)]);
grid on
hold on
errorbar(phs_obs, f_obs, phs_obs_err, 'horizontal', 'ko')
% legend('光滑反演正演响应', '观测数据')
% title('光滑反演正演响应')
xticks([0, 30, 60, 90])
yticks(logspace(log10(f_obs(1)), log10(f_obs(end)), log10(f_obs(end))*2+1))
set(gca,'FontName','Times New Roman', 'FontSize', 10)
xlabel('Phase/°')
ylabel('Frequency/Hz')
hold off