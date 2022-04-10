function lh = likelihood_func(Cd, rho_obs_log, phs_obs, f_obs, rho_log, z_log, k_weight, k_err)
%LIKELIHOOD_FUNC 似然函数
% 

[rho_log, phs] = forward_func(rho_log, z_log, f_obs);

rho_obs_log_std = std(rho_obs_log);
phs_obs_std = std(phs_obs);
rho_obs_log_mean = mean(rho_obs_log);
phs_obs_mean = mean(phs_obs);

% 标准化数据
rho_obs_log_standardize = (rho_obs_log-rho_obs_log_mean)/rho_obs_log_std;
phs_obs_standardize = (phs_obs-phs_obs_mean)/phs_obs_std;
rho_log_standardize = (rho_log-rho_obs_log_mean)/rho_obs_log_std;
phs_standardize = (phs-phs_obs_mean)/phs_obs_std;

phi = norm((Cd^(-1/2))*([rho_log_standardize; phs_standardize.*k_weight]-[rho_obs_log_standardize; phs_obs_standardize.*k_weight]) .* k_err)^2;
% phi = sum(((Cd^(-1/2))*([rho_log; phs.*k_phs]-[rho_obs_log; phs_obs.*k_phs]) .* k_err).^2);
% phi = norm((Cd^(-1/2))*([d_log; phs.*k_phs]-[d_obs_log; phs_obs.*k_phs]) .* k_err);
% phi = norm((Cd^(-1/2))*(d_log-d_obs_log));

lh = exp(-phi/2);
% lh = exp(-phi/2);

end

