function ovar = proposal_func(ilgv, flag, var_mesh)
%PROPOSAL_FUNC 建议分布函数
% var 需要扰动的参数 flag 扰动参数类型（rho、z）

olgv = -1;
while olgv <= 0
    % 均匀分布方法，rho~U(-1,1)，z~U(-0.5,0.5)
    if strcmp(flag, 'rho')
        olgv = ilgv + 2*rand-1;
    elseif strcmp(flag, 'z')
        olgv = ilgv + rand*0.5-0.25;
    end
    
    ovar = mesh_func(olgv, var_mesh);
end

end
