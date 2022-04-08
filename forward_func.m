function [rho_a_log, phase] = forward_func(m_log, z_log, f_obs)
%FORWARD_FUNC ��ص��һά���ݳ���
% m ģ�� z �½������ d ������Ӧ����
m = 10.^m_log;
z = 10.^z_log;

mu0 = 4*pi*10^-7; % ��մŵ���
h = [z(1); z(2:end)-z(1:end-1)]; % ����

R = ones(length(f_obs), 1);

for j = (length(m)-1):-1:1 % ����
    Ri = sqrt(m(j+1)/m(j))*R;
    L = (1-Ri)./(1+Ri);
    Li = L .* exp(-2*sqrt(-2i*pi*f_obs*mu0/m(j))*h(j));
    R = (1-Li)./(1+Li);
end
rho_a_log = log10(m(1)*abs(R).^2); % �ӵ�����
Z = R .* -sqrt(2i*pi*f_obs*mu0*m(1));
% phase = atan(imag(Z)./real(Z))*180/pi; % ��λ
phase = -atan2(real(Z), -imag(Z))*180/pi; % �˴��㷨������������ģ����ǽ����ȷ
end