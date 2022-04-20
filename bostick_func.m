function [rho_bostick, z_bostick] = bostick_func(rho_obs_log, f_obs)
%BOSTICK_FUNC bostick反演生成初始模型
% 
mu0 = 4*pi*10^-7; % 真空磁导率

z_bostick = log10(sqrt((10.^rho_obs_log)./(2*pi*f_obs*mu0)));
m_f = gradient(rho_obs_log) ./ log10(1./f_obs);
rho_bostick = real(rho_obs_log + log10((1+m_f)./(1-m_f)));
end

