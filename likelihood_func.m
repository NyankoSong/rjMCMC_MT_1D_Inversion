function [lh, L2_err] = likelihood_func(Cd, rhoa_obs_log, phs_obs, f_rhoa, f_phs, rho_log, z_log, k_weight, data_scale_factor)
%LIKELIHOOD_FUNC 似然函数
% 

[rhoa_log, ~] = forward_func(rho_log, z_log, f_rhoa);
[~, phs] = forward_func(rho_log, z_log, f_phs);

L2_err = norm((Cd^(-1/2))*([rhoa_log.*(2-k_weight); phs.*k_weight]-[rhoa_obs_log.*(2-k_weight); phs_obs.*k_weight]))^2;

% lh = (2*pi*data_scale_factor^2)^(-length(Cd)/2) * det(Cd)^(-1/2) * exp(-([rhoa_log.*(2-k_weight); phs.*k_weight]-[rhoa_obs_log.*(2-k_weight); phs_obs.*k_weight])' * (data_scale_factor^2*Cd)^(-1)*([rhoa_log.*(2-k_weight); phs.*k_weight]-[rhoa_obs_log.*(2-k_weight); phs_obs.*k_weight])/2);
lh = exp(-norm((data_scale_factor^2*Cd)^(-1/2)*([rhoa_log.*(2-k_weight); phs.*k_weight]-[rhoa_obs_log.*(2-k_weight); phs_obs.*k_weight]))^2/2);
end

