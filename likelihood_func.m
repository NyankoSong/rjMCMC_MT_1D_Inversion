function lh = likelihood_func(Cd, d_obs_log, phs_obs, f_obs, m_log, z_log, k_phs, k_err)
%LIKELIHOOD_FUNC ��Ȼ����
% lh ��Ȼ��

[d_log, phs] = forward_func(m_log, z_log, f_obs);

% TODO:������λȨ�ز���
phi = norm((Cd^(-1/2))*([d_log; phs.*k_phs]-[d_obs_log; phs_obs.*k_phs]) .* k_err);
% phi = norm((Cd^(-1/2))*(d_log-d_obs_log));

lh = exp(-phi/2);

end

