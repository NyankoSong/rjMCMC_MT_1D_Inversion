function [lh, phi_sqrt] = likelihood_func(Cd, rhoa_obs_log, phs_obs, f_rhoa, f_phs, rho_log, z_log, k_weight, scale_factor)
%LIKELIHOOD_FUNC 似然函数
% 

[rhoa_log, ~] = forward_func(rho_log, z_log, f_rhoa);
[~, phs] = forward_func(rho_log, z_log, f_phs);

% 标准化数据
% rhoa_obs_log_std = std(rhoa_obs_log);
% phs_obs_std = std(phs_obs);
% rhoa_obs_log_mean = mean(rhoa_obs_log);
% phs_obs_mean = mean(phs_obs);
% 
% rhoa_obs_log_standardize = (rhoa_obs_log-rhoa_obs_log_mean)/rhoa_obs_log_std;
% phs_obs_standardize = (phs_obs-phs_obs_mean)/phs_obs_std;
% rhoa_log_standardize = (rhoa_log-rhoa_obs_log_mean)/rhoa_obs_log_std;
% phs_standardize = (phs-phs_obs_mean)/phs_obs_std;

% phi = norm((Cd^(-1/2))*([rhoa_log_standardize.*(2-k_weight); phs_standardize.*k_weight]-[rhoa_obs_log_standardize.*(2-k_weight); phs_obs_standardize.*k_weight]) .* k_err_vec)^2;
% phi = norm((Cd^(-1/2))*([rhoa_log.*(2-k_weight); phs.*k_weight]-[rhoa_obs_log.*(2-k_weight); phs_obs.*k_weight]))^2 * k_err;
phi_sqrt = norm(((scale_factor^2*Cd)^(-1/2))*([rhoa_log.*(2-k_weight); phs.*k_weight]-[rhoa_obs_log.*(2-k_weight); phs_obs.*k_weight]));
phi = phi_sqrt^2;
% phi = sum(((Cd^(-1/2))*([rho_log.*(2-k_weight); phs.*k_weight]-[rho_obs_log.*(2-k_weight); phs_obs.*k_weight]) .* k_err_vec).^2);
% phi = norm((Cd^(-1/2))*([d_log; phs.*k_phs]-[d_obs_log; phs_obs.*k_phs]) .* k_err);
% phi = norm((Cd^(-1/2))*(d_log-d_obs_log));

lh = exp(-phi/2);
% lh = exp(-phi/2);

end

