function [rho_a_log, phase] = forward_func(m_log, z_log, f_obs)
%FORWARD_FUNC 大地电磁一维正演程序
% m 模型 z 下界面深度 d 正演响应向量
m = 10.^m_log;
z = 10.^z_log;

mu0 = 4*pi*10^-7; % 真空磁导率
h = [z(1); z(2:end)-z(1:end-1)]; % 层厚度

R = ones(length(f_obs), 1);

for j = (length(m)-1):-1:1 % 正演
    Ri = sqrt(m(j+1)/m(j))*R;
    L = (1-Ri)./(1+Ri);
    Li = L .* exp(-2*sqrt(-2i*pi*f_obs*mu0/m(j))*h(j));
    R = (1-Li)./(1+Li);
end
rho_a_log = log10(m(1)*abs(R).^2); % 视电阻率
Z = R .* -sqrt(2i*pi*f_obs*mu0*m(1));
% phase = atan(imag(Z)./real(Z))*180/pi; % 相位
phase = -atan2(real(Z), -imag(Z))*180/pi; % 此处算法可能是有问题的，但是结果正确
end