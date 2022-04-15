function ovar = perturb_func(ilgv, flag, var_mesh)
%PROPOSAL_FUNC 扰动函数
% ilgv 需要扰动的参数 flag 扰动参数类型（rho、z）

olgv = -1;
while olgv <= 0
    % 均匀分布方法，rho~U(-1,1)，z~U(-0.5,0.5)
    if strcmp(flag, 'rho')
        olgv = ilgv + 2*rand-1;
    elseif strcmp(flag, 'z')
        olgv = ilgv + rand*0.5-0.25;
    end
    
%     % 高斯分布方法，rho~N(0,1)，z~N(0,0.5)
%     if strcmp(flag, 'rho')
%         olgv = ilgv + randn;
%     elseif strcmp(flag, 'z')
%         olgv = ilgv + randn*0.5;
%     end
    
    ovar = mesh_func(olgv, var_mesh);
end

end
