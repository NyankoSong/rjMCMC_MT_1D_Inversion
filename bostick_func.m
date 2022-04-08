function [m_bostick, z_bostick] = bostick_func(rho_obs_log, phase_obs, f_obs)
%BOSTICK_FUNC bostick反演生成初始模型
% 
mu0 = 4*pi*10^-7; % 真空磁导率
rho_obs = 10.^rho_obs_log;

z_bostick = log10(sqrt(rho_obs./(2*pi*f_obs*mu0)));
m_bostick = log10(abs(rho_obs.*(180./(2*phase_obs)-1)));
end

