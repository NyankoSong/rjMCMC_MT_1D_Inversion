function ovar = proposal_func(ilgv, flag, var_mesh)
%PROPOSAL_FUNC ����ֲ�����
% var ��Ҫ�Ŷ��Ĳ��� flag �Ŷ��������ͣ�rho��z��

% �����Ŷ�
olgv = -1;
while olgv <= 0
    if strcmp(flag, 'rho')
        olgv = ilgv + 2*rand-1;
    elseif strcmp(flag, 'z')
        olgv = ilgv + rand*0.5-0.25;
    end
    
    ovar = mesh_func(olgv, var_mesh);
end

